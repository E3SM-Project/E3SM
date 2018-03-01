
      subroutine CLS_MAPS( )
!-----------------------------------------------------------------------
!        ... Form the individual method reaction maps from the
!            overall reaction map
!-----------------------------------------------------------------------

      use IO, only      : lout
      use VAR_MOD, only : extcnt
      use RXT_MOD, only : cls_rxt_map, cls_rxt_cnt, rxmap, rxmcnt, &
			  prdmap, prdcnt, hetmap, hetcnt, usrmap, &
			  usrcnt, prd_lim, prd_limp1
      use MO_MAKE_MAP, only : MAKE_MAP
      use MO_PRD_MAP,  only : PRD_MAP

      implicit none
     
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, allocatable  ::   template(:,:)
      integer  ::   i, j, astat
      integer  ::   class
      integer  ::   index
      integer  ::   row
      integer  ::   rxno
      integer  ::   rxtnt(2), rxtnt_cls(2)
      
      integer  ::  XLATE
      
!-----------------------------------------------------------------------
!   	... In the following the 1st column of template
!           represents the count of products in each class.
!           The 2nd column represents the class # of the product.
!           On input the 3rd column represents the master
!           product number and on output represents the product
!           number in the specific class.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        ... Scan the base "pure" production map
!            for class pure production entries
!-----------------------------------------------------------------------
      ALLOCATE( template(MAX(5,prd_lim),3),stat=astat )
      if( astat /= 0 ) then
         write(lout,*) 'CLS_MAP: Failed to allocate template array; error = ',astat
	 stop
      end if
      do i = 1,prdcnt
         template(:prd_lim,3) = prdmap(i,2:prd_limp1)
         call PRD_MAP( template )
	 rxno = prdmap(i,1)
         do class = 1,5
            if( template(class,1) /= 0 ) then
	       index = cls_rxt_cnt(1,class) + 1
               call MAKE_MAP( cls_rxt_map(index,:,class), &
                              cls_rxt_cnt(1,class), &
                              class, &
                              rxno, &
                              template(class,1), &
                              template )
            end if
         end do
      end do

!-----------------------------------------------------------------------
!        ... Scan the base linear reaction map
!            for class pure production entries
!-----------------------------------------------------------------------
      do i = 1,rxmcnt(1)
         template(:prd_lim,3) = rxmap(i,3:prd_lim+2,1)
         call PRD_MAP( template )
         rxtnt(1)     = ABS( rxmap(i,2,1) )
         rxtnt_cls(1) = XLATE( rxtnt(1) )
	 rxno = rxmap(i,1,1)
         do class = 1,5
            if( class /= rxtnt_cls(1) .and. template(class,1) /= 0 ) then
	       index = cls_rxt_cnt(1,class) + 1
               call MAKE_MAP( cls_rxt_map(index,:,class), &
                              cls_rxt_cnt(1,class), &
                              class, &
                              rxno, &
                              template(class,1), &
                              template )
               cls_rxt_map(index,2,class) = ABS( rxmap(i,2,1) )
            end if
         end do
      end do

!-----------------------------------------------------------------------
!        ... Scan the base nonlinear reaction map
!            for class pure production entries
!-----------------------------------------------------------------------
      do i = 1,rxmcnt(2)
         template(:prd_lim,3) = rxmap(i,4:prd_lim+3,2)
         call PRD_MAP( template )
         rxtnt     = ABS( rxmap(i,2:3,2) )
	 rxno = rxmap(i,1,2)
         do j = 1,2
            rxtnt_cls(j) = XLATE( rxtnt(j) )
         end do
         do class = 1,5
            if( class /= rxtnt_cls(1) .and. class /= rxtnt_cls(2) .and. template(class,1) /= 0 ) then
	       index = cls_rxt_cnt(1,class) + 1
               call MAKE_MAP( cls_rxt_map(index,:,class), &
                              cls_rxt_cnt(1,class), &
                              class, &
                              rxno, &
                              template(class,1), &
                              template )
               cls_rxt_map(index,2:3,class) = ABS( rxmap(i,2:3,2) )
            end if
         end do
      end do
      
!-----------------------------------------------------------------------
!        ... Scan the base linear reaction map
!            for entries in the class linear map
!-----------------------------------------------------------------------
      do i = 1,rxmcnt(1)
         template(:prd_lim,3) = rxmap(i,3:prd_lim+2,1)
         call PRD_MAP( template )
         rxtnt(1)     = ABS( rxmap(i,2,1) )
         class        = XLATE( rxtnt(1) )
	 rxno         = rxmap(i,1,1)
         if( template(class,1) /= 0 ) then
            index = MAX( SUM(cls_rxt_cnt(1:2,class))+1,1 )
            call MAKE_MAP( cls_rxt_map(index,:,class), &
                           cls_rxt_cnt(2,class), &
                           class, &
                           rxno, &
                           template(class,1), &
                           template )
            cls_rxt_map(index,2,class) = rxmap(i,2,1)
         else if( rxmap(i,2,1) > 0 ) then
	    cls_rxt_cnt(2,class)       = cls_rxt_cnt(2,class) + 1
            row                        = SUM( cls_rxt_cnt(1:2,class) )
            cls_rxt_map(row,1:2,class) = rxmap(i,1:2,1)
         end if 
      end do
      
!-----------------------------------------------------------------------
!        ... Scan the base nonlinear reaction map
!            for entries in the class linear map
!-----------------------------------------------------------------------
      do i = 1,rxmcnt(2)
         do j = 1,2
            rxtnt(j)     = ABS( rxmap(i,j+1,2) )
            rxtnt_cls(j) = XLATE( rxtnt(j) )
         end do
         if( rxtnt_cls(1) /= rxtnt_cls(2) ) then
            template(:prd_lim,3) = rxmap(i,4:prd_lim+3,2)
            call PRD_MAP( template )
            do j = 1,2
               class = rxtnt_cls(j)
               if( template(class,1) /= 0 ) then
                  index = MAX( SUM(cls_rxt_cnt(1:2,class))+1,1 )
                  call MAKE_MAP( cls_rxt_map(index,:,class), &
                                 cls_rxt_cnt(2,class), &
                                 class, &
                                 rxmap(i,1,2), &
                                 template(class,1), &
                                 template )
                  cls_rxt_map(index,2,class) = rxmap(i,j+1,2)
                  if( j == 1 ) then
                     cls_rxt_map(index,3,class) = ABS( rxmap(i,3,2) )
                  else
                     cls_rxt_map(index,3,class) = ABS( rxmap(i,2,2) )
                  end if
               else if( rxmap(i,j+1,2) > 0 ) then
		  cls_rxt_cnt(2,class) = cls_rxt_cnt(2,class) + 1
                  row = cls_rxt_cnt(1,class) + cls_rxt_cnt(2,class)
                  cls_rxt_map(row,1,class) = rxmap(i,1,2)
                  if( j == 1 ) then
                     cls_rxt_map(row,2:3,class) = ABS( rxmap(i,2:3,2) )
                  else
                     cls_rxt_map(row,2:3,class) = ABS( rxmap(i,3:2:-1,2) )
                  end if
               end if
            end do
         end if
      end do

!-----------------------------------------------------------------------
!        ... Scan the base nonlinear reaction map
!            for entries in the class nonlinear map
!-----------------------------------------------------------------------
      do i = 1,rxmcnt(2)
         do j = 1,2
            rxtnt(j)     = ABS( rxmap(i,j+1,2) )
            rxtnt_cls(j) = XLATE( rxtnt(j) )
         end do
         if( rxtnt_cls(1) == rxtnt_cls(2) ) then
            template(:prd_lim,3) = rxmap(i,4:prd_lim+3,2)
            call PRD_MAP( template )
            class = rxtnt_cls(1)
            if( template(class,1) /= 0 ) then
               index = MAX( 1,SUM(cls_rxt_cnt(1:3,class))+1 )
               call MAKE_MAP( cls_rxt_map(index,:,class), &
                              cls_rxt_cnt(3,class), &
                              class, &
                              rxmap(i,1,2), &
                              template(class,1), &
                              template )
               cls_rxt_map(index,2:3,class) = rxmap(i,2:3,2)
            else if( rxmap(i,2,2) > 0 .or. rxmap(i,3,2) > 0 ) then
	       cls_rxt_cnt(3,class) = cls_rxt_cnt(3,class) + 1
               row = SUM( cls_rxt_cnt(1:3,class) )
               cls_rxt_map(row,1:3,class) = rxmap(i,1:3,2)
            end if
         end if
      end do

!-----------------------------------------------------------------------
!        ... Scan the heterogeneous reactions
!-----------------------------------------------------------------------
      do i = 1,hetcnt
	 rxtnt(1)     = ABS( hetmap(i,1) )
         class        = XLATE( rxtnt(1) )
	 cls_rxt_cnt(4,class) = cls_rxt_cnt(4,class) + 1
         index = SUM( cls_rxt_cnt(1:4,class) )
	 cls_rxt_map(index,1,class) = i
	 cls_rxt_map(index,2,class) = rxtnt(1)
      end do

!-----------------------------------------------------------------------
!        ... Scan the extraneous forcing
!-----------------------------------------------------------------------
      do i = 1,usrcnt
	 rxtnt(1)      = ABS( usrmap(i) )
         class         = XLATE( rxtnt(1) )
	 extcnt(class) = extcnt(class) + 1
         index = SUM( cls_rxt_cnt(1:4,class) ) + extcnt(class)
	 cls_rxt_map(index,1,class) = i
	 cls_rxt_map(index,2,class) = rxtnt(1)
      end do

      DEALLOCATE( template )

      end subroutine CLS_MAPS
