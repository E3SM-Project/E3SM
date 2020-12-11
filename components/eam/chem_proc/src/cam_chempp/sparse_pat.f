
      subroutine SPARSITY_PAT( clscnt, &
                               clsmap, &
                               cls_rxt_cnt, &
                               cls_rxt_map, &
                               sparse_pat )
!-----------------------------------------------------------------------
!	... Set the jacobian matrix sparsity pattern
!-----------------------------------------------------------------------
     
      use VAR_MOD, only : var_lim
      use RXT_MOD, only : rxt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!
!            The columns of the cls_rxt_cnt represent the reaction count
!	     for each class with the following row conontation:
!		(1) - independent reactions
!		(2) - linear reactions
!		(3) - nonlinear reactions
!		(4) - heterogeneous processes
!-----------------------------------------------------------------------
      integer, intent(in) ::  clscnt, &
                              clsmap(var_lim), &
                              cls_rxt_map(rxt_lim,prd_lim+3), &
                              cls_rxt_cnt(4)
      logical, intent(out)::  sparse_pat(clscnt,clscnt)
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer  ::   i, k, kl, ku, l, m
      integer  ::   target
      integer  ::   species
      integer, allocatable  ::   indexer(:)
      logical, allocatable :: match_mask(:,:)
      logical, allocatable :: pmask(:,:)
      
      if( ALLOCATED( match_mask ) ) then
         DEALLOCATE( match_mask )
      end if
      if( ALLOCATED( pmask ) ) then
	 DEALLOCATE( pmask )
      end if
      if( ALLOCATED( indexer ) ) then
	 DEALLOCATE( indexer )
      end if
      k = SUM( cls_rxt_cnt(:) )
      ALLOCATE( match_mask(k,3) )
      ALLOCATE( indexer(k) )
      if( SUM( cls_rxt_cnt(2:3) ) /= 0 ) then
         ALLOCATE( pmask(k,prd_lim) )
      end if
      sparse_pat = .false.
      do i = 1,clscnt
         sparse_pat(i,i) = .true.                   ! assume only diagonal entries
      end do

Species_loop : &
      do species = 1,clscnt
!-----------------------------------------------------------------------
!   	... Check for non-linear losses
!-----------------------------------------------------------------------
         target = clsmap(species)
	 kl = SUM( cls_rxt_cnt(:2) ) + 1
	 ku = SUM( cls_rxt_cnt(:3) )
	 do i = 1,2
	    match_mask(kl:ku,i) = cls_rxt_map(kl:ku,i+1) == target
	    where( match_mask(kl:ku,i) )
	       indexer(kl:ku) = 6/(i+1)
	    endwhere
	 end do
         match_mask(kl:ku,1) = match_mask(kl:ku,1) .or. match_mask(kl:ku,2)
	 if( COUNT( match_mask(kl:ku,1) ) /= 0 ) then
            do k = kl,ku
               if( match_mask(k,1) ) then
	          m = ABS( cls_rxt_map(k,indexer(k)) )
		  if( m /= target ) then
		     do i = 1,clscnt
                        if( clsmap(i) == m ) then
			   sparse_pat(species,i) = .true.
		           exit
	                end if
		     end do
		  end if
	       end if
	    end do
         end if
!-----------------------------------------------------------------------
!   	... Check for production from linear and nonlinear reactions
!-----------------------------------------------------------------------
	 kl = cls_rxt_cnt(1) + 1
	 do k = kl,ku
            pmask(k,:) = cls_rxt_map(k,4:prd_lim+3) == species
	    match_mask(k,1) = ANY( pmask(k,:) )
	 end do
	 if( COUNT( match_mask(kl:ku,1) ) /= 0 ) then
	    do k = kl,ku
               if( match_mask(k,1) ) then
	          do i = 2,3
	             m = ABS( cls_rxt_map(k,i) )
		     if( m /= 0 ) then
		        do l = 1,clscnt
                           if( clsmap(l) == m ) then
			      sparse_pat(species,l) = .true.
			      exit
	                   end if
                        end do
                     end if
	          end do
	       end if
	    end do
	 end if
      end do Species_loop

      if( ALLOCATED( match_mask ) ) then
	 DEALLOCATE( match_mask )
      end if
      if( ALLOCATED( pmask ) ) then
	 DEALLOCATE( pmask )
      end if
      if( ALLOCATED( indexer ) ) then
	 DEALLOCATE( indexer )
      end if

      end subroutine SPARSITY_PAT
