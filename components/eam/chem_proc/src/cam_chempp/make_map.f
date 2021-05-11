
      module MO_MAKE_MAP

      CONTAINS

      subroutine MAKE_MAP( cls_rxt_map, &
                           cls_rxt_cnt, &
                           clsno, &
                           rxno, &
                           cls_prd_cnt, &
                           template )
     
      use RXT_MOD, only : rxt_lim, prd_lim

      implicit none

!------------------------------------------------------------------------
!	... Dummy args
!------------------------------------------------------------------------
      integer, intent(in)    ::   clsno,  rxno,  cls_prd_cnt
      integer, intent(in)    ::   template(:,:)
      integer, intent(inout) ::   cls_rxt_cnt
      integer, intent(inout) ::   cls_rxt_map(:)
     
!------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------
      integer  ::  count
      integer  ::  k, kp3
     
      count          = 0
      cls_rxt_cnt    = cls_rxt_cnt + 1
      cls_rxt_map(1) = rxno
      do k = 1,prd_lim
         kp3 = k + 3
         if( template(k,2) == clsno ) then
            count = count + 1
            cls_rxt_map(kp3) = template(k,3)
            if( count == cls_prd_cnt ) then
               exit
            end if
         else
            cls_rxt_map(kp3) = -99
         end if
      end do

      end subroutine MAKE_MAP

      end module MO_MAKE_MAP
