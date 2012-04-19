#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


  interface

#if 0
     subroutine initgridvertex(GridVertex)
       use gridgraph_mod
       use coordinate_systems_mod
       use cube_mod
       implicit none
       type (GridVertex_t), intent(out)    :: GridVertex(:)
     end subroutine initgridvertex

     subroutine initgridedge(GridEdge,GridVertex)
       use gridgraph_mod
       use schedule_mod
       implicit none
       type (GridEdge_t), intent(inout) :: GridEdge(:)
       type (GridVertex_t), intent(in)   :: GridVertex(:)
     end subroutine initgridedge

     subroutine initgridgraph(GridEdge,GridVertex)
       use gridgraph_mod
       use schedule_mod
       implicit none
       type (GridEdge_t), intent(out)     :: GridEdge(:)
       type (GridVertex_t), intent(out)   :: GridVertex(:)
     end subroutine initgridgraph

!     subroutine genchecksum(TestField,Checksum,GridEdge)
!       use gridgraph_mod
!       use field_mod
!       use control_mod
!       type (GridEdge_t), intent(in),target   :: GridEdge(:)
!       type (IntField_t), intent(inout),target  :: TestField(:)
!       type (IntField_t), intent(inout),target  :: Checksum(:)
!     end subroutine genchecksum

     subroutine genmetispart(GridEdge,GridVertex)
       use gridgraph_mod
       use schedule_mod
       use control_mod
       implicit none
       type (GridVertex_t), intent(inout) :: GridVertex(:)
       type (GridEdge_t),intent(inout)    :: GridEdge(:)
     end subroutine genmetispart

     subroutine genmetasched(MetaEdge,MetaVertex)
       use metagraph_mod
       use schedule_mod
!JMD       use output_mod
       implicit none
       type (MetaEdge_t),   intent(inout),target :: MetaEdge(:)
       type (MetaVertex_t), intent(inout),target :: MetaVertex(:)
     end subroutine genmetasched

     subroutine prmetasched()
       use schedule_mod
       implicit none
     end subroutine prmetasched

     subroutine testmetasched_nonblock() 
!JMD       use field_mod
       use metagraph_mod
       use schedule_mod
!JMD       use output_mod
       use control_mod
       use parallel_mod
       implicit none
     end subroutine testmetasched_nonblock

     subroutine testchecksum(par,GridEdge)
       use gridgraph_mod
       use field_mod
       use metagraph_mod
       use schedule_mod
       use control_mod
       use parallel_mod
       implicit none
       type (Parallel_t),     intent(in)         :: par
       type (GridEdge_t),     intent(in),target  :: GridEdge(:)
     end subroutine testchecksum
   
     subroutine readnamelist()
       use dimensions_mod
       use parallel_mod
       use control_mod
       implicit none
     end subroutine readnamelist
     
     subroutine  initMetaGraph(MetaVertex,MetaEdge,GridVertex,GridEdge)
       use gridgraph_mod
       use metagraph_mod
       use schedule_mod
       implicit none
       type (MetaVertex_t), intent(out),target :: MetaVertex(:)
       type (MetaEdge_t),   intent(out),target :: MetaEdge(:)
       type (GridVertex_t), intent(in),target  :: GridVertex(:)
       type (GridEdge_t),   intent(in),target  :: GridEdge(:)
     end subroutine initMetaGraph
#endif

   
  end interface
