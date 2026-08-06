#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  MT 2018/5 updated to use pack/unpack with variable nlyr 
!
module edge_mod
  use edge_mod_base, only: initLongEdgeBuffer, FreeLongEdgeBuffer, LongEdgeVpack, LongEdgeVunpackMIN, initEdgeBuffer, initEdgeSBuffer, FreeEdgeBuffer, &
                           edgeVpack, edgeVunpack, edgeVpack_nlyr, edgeVunpack_nlyr,       &
                           edgeVunpackMIN, edgeVunpackMAX, edgeDGVpack, edgeDGVunpack, edgeVunpackVert, edgeDefaultVal, initGhostBuffer3D, FreeGhostBuffer3D, &
                           ghostVpackfull, ghostVunpackfull, ghostVpack_unoriented, ghostVunpack_unoriented, ghostVpack3d, ghostVunpack3d, &
                           edgeSpack, edgeSunpackMin, edgeSunpackMax, edge_g
  use kinds, only : int_kind, log_kind, real_kind
  use dimensions_mod, only : max_neigh_edges, nelemd, np
  use perf_mod, only: t_startf, t_stopf, t_adj_detailf ! _EXTERNAL
  use coordinate_systems_mod, only : cartesian3D_t
  use schedtype_mod, only : cycle_t, schedule_t, schedule
  use parallel_mod, only : abortmp, haltmp, MPIreal_t, iam,parallel_t, &
      MAX_ACTIVE_MSG, HME_status_size, BNDRY_TAG_BASE
  use edgetype_mod, only : edgedescriptor_t, edgebuffer_t, &
      Longedgebuffer_t, initedgebuffer_callid
  use element_mod, only : element_t
  implicit none
  private

  public :: initLongEdgeBuffer, FreeLongEdgeBuffer, LongEdgeVpack, LongEdgeVunpackMIN, initEdgeBuffer, initEdgeSBuffer, FreeEdgeBuffer,&
       edgeVpack, edgeVunpack, edgeVpack_nlyr, edgeVunpack_nlyr,      &
       edgeVunpackMIN, edgeVunpackMAX, edgeDGVpack, edgeDGVunpack, edgeVunpackVert, edgeDefaultVal, initGhostBuffer3D, FreeGhostBuffer3D, &
       ghostVpackfull, ghostVunpackfull, ghostVpack_unoriented, ghostVunpack_unoriented, ghostVpack3d, ghostVunpack3d, &
       edgeSpack, edgeSunpackMin, edgeSunpackMax, edge_g
  public :: edgeSpack_openacc
  public :: edgeSunpackMin_openacc
  public :: edgeSunpackMax_openacc
  public :: edgeVpack_openacc
  public :: edgeVunpack_openacc


contains

  subroutine edgeSpack_openacc(edge,v,vlyr,kptr,nlyr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : max_corner_elem
    use control_mod   , only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod      , only : t_startf, t_stopf
    use parallel_mod  , only : haltmp
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    type(EdgeBuffer_t)     ,intent(inout) :: edge
    integer                ,intent(in   ) :: vlyr,nlyr
    integer                ,intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    real (kind=real_kind)  ,intent(in   ) :: v(vlyr,tdim,nelemd)
    ! Local variables
    type (EdgeDescriptor_t),pointer            :: desc        ! =>elem(ie)%desc
    integer :: i,k,ir,ll,is,ie,in,iw,el,kc,kk
    integer, parameter :: kchunk = 64
    call t_startf('edge_s_pack')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeSpack: Buffer overflow1: size of the vertical dimension must be increased!')
    if (edge%nlyr_max < (kptr+vlyr) ) call haltmp('edgeSpack: Buffer overflow: size of the vertical dimension must be increased!')
    edge%nlyr=nlyr      ! set total amount packed, for use by bndry_exchange

    !$acc parallel loop gang collapse(2) present(v,edge) vector_length(kchunk)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) then
            edge%buf(nlyr*edge%desc(el)%putmapS(south)+kptr+k) = v(k,tl,el)
            edge%buf(nlyr*edge%desc(el)%putmapS(east )+kptr+k) = v(k,tl,el)
            edge%buf(nlyr*edge%desc(el)%putmapS(north)+kptr+k) = v(k,tl,el)
            edge%buf(nlyr*edge%desc(el)%putmapS(west )+kptr+k) = v(k,tl,el)
          endif
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(edge%desc(el)%putmapS(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapS(ll)+max_corner_elem*(kptr+k-1)+i) = v(k,tl,el)
              ll = swest+1*max_corner_elem+i-1; if(edge%desc(el)%putmapS(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapS(ll)+max_corner_elem*(kptr+k-1)+i) = v(k,tl,el)
              ll = swest+2*max_corner_elem+i-1; if(edge%desc(el)%putmapS(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapS(ll)+max_corner_elem*(kptr+k-1)+i) = v(k,tl,el)
              ll = swest+3*max_corner_elem+i-1; if(edge%desc(el)%putmapS(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapS(ll)+max_corner_elem*(kptr+k-1)+i) = v(k,tl,el)
            endif
          enddo
        enddo
      enddo
    enddo
    call t_stopf('edge_s_pack')
  end subroutine edgeSpack_openacc

  subroutine edgeSunpackMin_openacc(edge,v,vlyr,kptr,nlyr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    type(EdgeBuffer_t)    , intent(in   ) :: edge
    integer               , intent(in   ) :: vlyr,nlyr
    integer               , intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    real(kind=real_kind)  , intent(inout) :: v(vlyr,tdim,nelemd)
    ! Local
    type (EdgeDescriptor_t),pointer            :: desc        ! =>elem(ie)%desc
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk
    integer, parameter :: kchunk = 64
    real(kind=real_kind) :: vtmp(kchunk)
    call t_startf('edge_s_unpack_min')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeSunpackMin: Buffer overflow1: size of the vertical dimension must be increased!')
    if (edge%nlyr_max < (kptr+vlyr) ) call haltmp('edgeSunpackMin: Buffer overflow: size of the vertical dimension must be increased!')
    !$acc parallel loop gang collapse(2) present(v,edge) private(vtmp) vector_length(kchunk)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc cache(vtmp)
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k > vlyr) k = vlyr
          vtmp(kk) = v(k,tl,el)
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) then
            vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(south)+kptr+k) )
            vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(east )+kptr+k) )
            vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(north)+kptr+k) )
            vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(west )+kptr+k) )
          endif
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+1*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+2*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+3*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = min( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
            endif
          enddo
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) v(k,tl,el) = vtmp(kk)
        enddo
      enddo
    enddo
    call t_stopf('edge_s_unpack_min')
  end subroutine edgeSunpackMin_openacc

  subroutine edgeSunpackMax_openacc(edge,v,vlyr,kptr,nlyr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    type(EdgeBuffer_t)    , intent(in   ) :: edge
    integer               , intent(in   ) :: vlyr,nlyr
    integer               , intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    real(kind=real_kind)  , intent(inout) :: v(vlyr,tdim,nelemd)
    ! Local
    type (EdgeDescriptor_t),pointer            :: desc        ! =>elem(ie)%desc
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk
    integer, parameter :: kchunk = 64
    real(kind=real_kind) :: vtmp(kchunk)
    call t_startf('edge_s_unpack_max')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeSunpackMax: Buffer overflow1: size of the vertical dimension must be increased!')
    if (edge%nlyr_max < (kptr+vlyr) ) call haltmp('edgeSunpackMax: Buffer overflow: size of the vertical dimension must be increased!')
    !$acc parallel loop gang collapse(2) present(v,edge) private(vtmp) vector_length(kchunk)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc cache(vtmp)
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k > vlyr) k = vlyr
          vtmp(kk) = v(k,tl,el)
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) then
            vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(south)+kptr+k) )
            vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(east )+kptr+k) )
            vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(north)+kptr+k) )
            vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(west )+kptr+k) )
          endif
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+1*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+2*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
              ll = swest+3*max_corner_elem+i-1; if(edge%desc(el)%getmapS(ll) /= -1) vtmp(kk) = max( vtmp(kk) , edge%receive(nlyr*edge%desc(el)%getmapS(ll)+max_corner_elem*(kptr+k-1)+i) )
            endif
          enddo
        enddo
        !$acc loop vector
        do kk = 1 , kchunk
          k = (kc-1)*kchunk+kk
          if (k <= vlyr) v(k,tl,el) = vtmp(kk)
        enddo
      enddo
    enddo
    call t_stopf('edge_s_unpack_max')
  end subroutine edgeSunpackMax_openacc

  subroutine edgeVpack_openacc(edge,v,vlyr,kptr,nlyr,nets,nete,tdim,tl)
    use dimensions_mod, only : max_corner_elem
    use control_mod   , only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod      , only : t_startf, t_stopf
    use parallel_mod  , only : haltmp
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    type(EdgeBuffer_t)     ,intent(inout) :: edge
    integer                ,intent(in   ) :: vlyr,nlyr
    integer                ,intent(in   ) :: kptr
!    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    real (kind=real_kind)  ,intent(in   ) :: v(np,np,vlyr,tdim,nelemd)
    ! Local variables
    type (EdgeDescriptor_t),pointer            :: desc        ! =>elem(ie)%desc
    integer :: i,k,ir,ll,is,ie,in,iw,el,kc,kk
    integer, parameter :: kchunk = 32
    call t_startf('edge_pack')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeVpack: Buffer overflow1: size of the vertical dimension must be increased!')
    if (edge%nlyr_max < nlyr ) call haltmp('edgeVpack: Buffer overflow2: size of the vertical dimension must be increased!')
    edge%nlyr=nlyr  ! set total amount packed for use by bndry_exchange
    !$acc parallel loop gang collapse(2) present(v,edge) vector_length(kchunk*np)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , np
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              edge%buf(nlyr*edge%desc(el)%putmapP(south)+np*(kptr+k-1)+i) = v(i ,1 ,k,tl,el)
              edge%buf(nlyr*edge%desc(el)%putmapP(east )+np*(kptr+k-1)+i) = v(np,i ,k,tl,el)
              edge%buf(nlyr*edge%desc(el)%putmapP(north)+np*(kptr+k-1)+i) = v(i ,np,k,tl,el)
              edge%buf(nlyr*edge%desc(el)%putmapP(west )+np*(kptr+k-1)+i) = v(1 ,i ,k,tl,el)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , np
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ir = np-i+1
              if(edge%desc(el)%reverse(south)) edge%buf(nlyr*edge%desc(el)%putmapP(south)+np*(kptr+k-1)+ir) = v(i ,1 ,k,tl,el)
              if(edge%desc(el)%reverse(east )) edge%buf(nlyr*edge%desc(el)%putmapP(east )+np*(kptr+k-1)+ir) = v(np,i ,k,tl,el)
              if(edge%desc(el)%reverse(north)) edge%buf(nlyr*edge%desc(el)%putmapP(north)+np*(kptr+k-1)+ir) = v(i ,np,k,tl,el)
              if(edge%desc(el)%reverse(west )) edge%buf(nlyr*edge%desc(el)%putmapP(west )+np*(kptr+k-1)+ir) = v(1 ,i ,k,tl,el)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(edge%desc(el)%putmapP(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapP(ll)+max_corner_elem*(kptr+k-1)+i) = v(1 ,1 ,k,tl,el)
              ll = swest+1*max_corner_elem+i-1; if(edge%desc(el)%putmapP(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapP(ll)+max_corner_elem*(kptr+k-1)+i) = v(np,1 ,k,tl,el)
              ll = swest+2*max_corner_elem+i-1; if(edge%desc(el)%putmapP(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapP(ll)+max_corner_elem*(kptr+k-1)+i) = v(1 ,np,k,tl,el)
              ll = swest+3*max_corner_elem+i-1; if(edge%desc(el)%putmapP(ll) /= -1) edge%buf(nlyr*edge%desc(el)%putmapP(ll)+max_corner_elem*(kptr+k-1)+i) = v(np,np,k,tl,el)
            endif
          enddo
        enddo
      enddo
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack_openacc

  subroutine edgeVunpack_openacc(edge,v,vlyr,kptr,nlyr,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    type(EdgeBuffer_t)    , intent(in   ) :: edge
    integer               , intent(in   ) :: vlyr
    integer               , intent(in   ) :: kptr, nlyr
    !type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    ! Local
    type (EdgeDescriptor_t),pointer            :: desc        ! =>elem(ie)%desc
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk,glob_k,loc_ind,ii,jj, j
    integer, parameter :: kchunk = 32
    real(kind=real_kind) :: vtmp(np,np,kchunk)
    call t_startf('edge_unpack')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeVunpack: Buffer overflow1: size of the vertical dimension must be increased!')
    if (edge%nlyr_max < nlyr ) call haltmp('edgeVunpack: Buffer overflow2: size of the vertical dimension must be increased!')
    !$acc parallel loop gang collapse(2) present(v,edge) private(vtmp)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc cache(vtmp)
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              loc_ind = ((j-1)*np+i-1)*kchunk+kk-1
              ii = modulo(loc_ind,np)+1
              jj = modulo(loc_ind/np,np)+1
              k = loc_ind/np/np+1
              glob_k = (kc-1)*kchunk+k
              if (glob_k > vlyr) glob_k = vlyr
              vtmp(ii,jj,k) = v(ii,jj,glob_k,tl,el)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , np
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              vtmp(i ,1 ,kk) = vtmp(i ,1 ,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(south)+np*(kptr+k-1)+i)
              vtmp(np,i ,kk) = vtmp(np,i ,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(east )+np*(kptr+k-1)+i)
              vtmp(i ,np,kk) = vtmp(i ,np,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(north)+np*(kptr+k-1)+i)
              vtmp(1 ,i ,kk) = vtmp(1 ,i ,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(west )+np*(kptr+k-1)+i)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(edge%desc(el)%getmapP(ll) /= -1) vtmp(1  ,1 ,kk) = vtmp(1 ,1 ,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(ll)+max_corner_elem*(kptr+k-1)+i)
              ll = swest+1*max_corner_elem+i-1; if(edge%desc(el)%getmapP(ll) /= -1) vtmp(np ,1 ,kk) = vtmp(np,1 ,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(ll)+max_corner_elem*(kptr+k-1)+i)
              ll = swest+2*max_corner_elem+i-1; if(edge%desc(el)%getmapP(ll) /= -1) vtmp(1  ,np,kk) = vtmp(1 ,np,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(ll)+max_corner_elem*(kptr+k-1)+i)
              ll = swest+3*max_corner_elem+i-1; if(edge%desc(el)%getmapP(ll) /= -1) vtmp(np ,np,kk) = vtmp(np,np,kk) + edge%receive(nlyr*edge%desc(el)%getmapP(ll)+max_corner_elem*(kptr+k-1)+i)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              loc_ind = ((j-1)*np+i-1)*kchunk+kk-1
              ii = modulo(loc_ind,np)+1
              jj = modulo(loc_ind/np,np)+1
              k = loc_ind/np/np+1
              glob_k = (kc-1)*kchunk+k
              if (glob_k <= vlyr) v(ii,jj,glob_k,tl,el) = vtmp(ii,jj,k)
            enddo
          enddo
        enddo
      enddo
    enddo
    call t_stopf('edge_unpack')
  end subroutine edgeVunpack_openacc


end module edge_mod

