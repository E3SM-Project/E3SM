
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module edge_openacc_mod
  use kinds, only: real_kind, int_kind, log_kind
  use dimensions_mod, only: max_neigh_edges,nelemd,np,max_corner_elem
  implicit none
  private
  integer(kind=int_kind), allocatable :: putmapP(:,:)
  integer(kind=int_kind), allocatable :: getmapP(:,:)
  logical(kind=log_kind), allocatable :: reverse(:,:)
  integer :: nbuf
  logical :: maps_allocated = .false.

  public :: edgeVpack
  public :: edgeVunpack
  public :: edgeVunpackMin
  public :: edgeVunpackMax

contains

  subroutine alloc_maps(elem)
    use element_mod, only: element_t
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer :: ie
    nbuf = 4*(np+max_corner_elem)*nelemd

    allocate(putmapP(max_neigh_edges,nelemd))
    allocate(getmapP(max_neigh_edges,nelemd))
    allocate(reverse(max_neigh_edges,nelemd))
    do ie = 1 , nelemd
      putmapP(:,ie) = elem(ie)%desc%putmapP
      getmapP(:,ie) = elem(ie)%desc%getmapP
      reverse(:,ie) = elem(ie)%desc%reverse
    enddo
    !$acc enter data pcreate(putmapP,getmapP,reverse)
    !$acc update device(putmapP,getmapP,reverse)

    maps_allocated = .true.
  end subroutine alloc_maps

  subroutine edgeVpack(edgebuf,nlyr,v,vlyr,kptr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : max_corner_elem
    use control_mod   , only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod      , only : t_startf, t_stopf
    use parallel_mod  , only : haltmp
    use element_mod   , only : Element_t
    use edgetype_mod  , only : EdgeBuffer_t
    real(kind=real_kind)  , intent(inout) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer                ,intent(in   ) :: vlyr
    real (kind=real_kind)  ,intent(in   ) :: v(np,np,vlyr,tdim,nelemd)
    integer                ,intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    ! Local variables
    integer :: i,k,ir,ll,is,ie,in,iw,el,kc,kk
    integer, parameter :: kchunk = 32
    if (.not. maps_allocated) call alloc_maps(elem)
    call t_startf('edge_pack')
    if (nlyr < (kptr+vlyr) ) call haltmp('edgeVpack: Buffer overflow: size of the vertical dimension must be increased!')
    !$acc parallel loop gang collapse(2) present(v,putmapP,reverse,edgebuf) vector_length(kchunk*np)
    do el = nets , nete
      do kc = 1 , vlyr/kchunk+1
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , np
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              edgebuf(kptr+k,putmapP(south,el)+i) = v(i ,1 ,k,tl,el)
              edgebuf(kptr+k,putmapP(east ,el)+i) = v(np,i ,k,tl,el)
              edgebuf(kptr+k,putmapP(north,el)+i) = v(i ,np,k,tl,el)
              edgebuf(kptr+k,putmapP(west ,el)+i) = v(1 ,i ,k,tl,el)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , np
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ir = np-i+1
              if(reverse(south,el)) edgebuf(kptr+k,putmapP(south,el)+ir) = v(i ,1 ,k,tl,el)
              if(reverse(east ,el)) edgebuf(kptr+k,putmapP(east ,el)+ir) = v(np,i ,k,tl,el)
              if(reverse(north,el)) edgebuf(kptr+k,putmapP(north,el)+ir) = v(i ,np,k,tl,el)
              if(reverse(west ,el)) edgebuf(kptr+k,putmapP(west ,el)+ir) = v(1 ,i ,k,tl,el)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edgebuf(kptr+k,putmapP(ll,el)+1) = v(1 ,1 ,k,tl,el)
              ll = swest+1*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edgebuf(kptr+k,putmapP(ll,el)+1) = v(np,1 ,k,tl,el)
              ll = swest+2*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edgebuf(kptr+k,putmapP(ll,el)+1) = v(1 ,np,k,tl,el)
              ll = swest+3*max_corner_elem+i-1; if(putmapP(ll,el) /= -1) edgebuf(kptr+k,putmapP(ll,el)+1) = v(np,np,k,tl,el)
            endif
          enddo
        enddo
      enddo
    enddo
    call t_stopf('edge_pack')
  end subroutine edgeVpack

  subroutine edgeVunpack(edgebuf,nlyr,v,vlyr,kptr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk,glob_k,loc_ind,ii,jj, j
    integer, parameter :: kchunk = 32
    real(kind=real_kind) :: vtmp(np,np,kchunk)
    call t_startf('edge_unpack')
    !$acc parallel loop gang collapse(2) present(v,getmapP,edgebuf) private(vtmp)
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
              vtmp(i ,1 ,kk) = vtmp(i ,1 ,kk) + edgebuf(kptr+k,getmapP(south,el)+i)
              vtmp(np,i ,kk) = vtmp(np,i ,kk) + edgebuf(kptr+k,getmapP(east ,el)+i)
              vtmp(i ,np,kk) = vtmp(i ,np,kk) + edgebuf(kptr+k,getmapP(north,el)+i)
              vtmp(1 ,i ,kk) = vtmp(1 ,i ,kk) + edgebuf(kptr+k,getmapP(west ,el)+i)
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,1 ,kk) = vtmp(1 ,1 ,kk) + edgebuf(kptr+k,getmapP(ll,el)+1)
              ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,1 ,kk) = vtmp(np,1 ,kk) + edgebuf(kptr+k,getmapP(ll,el)+1)
              ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,np,kk) = vtmp(1 ,np,kk) + edgebuf(kptr+k,getmapP(ll,el)+1)
              ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,np,kk) = vtmp(np,np,kk) + edgebuf(kptr+k,getmapP(ll,el)+1)
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
  end subroutine edgeVunpack

  subroutine edgeVunpackMin(edgebuf,nlyr,v,vlyr,kptr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk,glob_k,loc_ind,ii,jj, j
    integer, parameter :: kchunk = 32
    real(kind=real_kind) :: vtmp(np,np,kchunk)
    call t_startf('edge_unpack_min')
    !$acc parallel loop gang collapse(2) present(v,getmapP,edgebuf) private(vtmp)
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
              vtmp(i ,1 ,kk) = min( vtmp(i ,1 ,kk) , edgebuf(kptr+k,getmapP(south,el)+i) )
              vtmp(np,i ,kk) = min( vtmp(np,i ,kk) , edgebuf(kptr+k,getmapP(east ,el)+i) )
              vtmp(i ,np,kk) = min( vtmp(i ,np,kk) , edgebuf(kptr+k,getmapP(north,el)+i) )
              vtmp(1 ,i ,kk) = min( vtmp(1 ,i ,kk) , edgebuf(kptr+k,getmapP(west ,el)+i) )
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,1 ,kk) = min( vtmp(1 ,1 ,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,1 ,kk) = min( vtmp(np,1 ,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,np,kk) = min( vtmp(1 ,np,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,np,kk) = min( vtmp(np,np,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
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
    call t_stopf('edge_unpack_min')
  end subroutine edgeVunpackMin

  subroutine edgeVunpackMax(edgebuf,nlyr,v,vlyr,kptr,elem,nets,nete,tdim,tl)
    use dimensions_mod, only : np, max_corner_elem
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use perf_mod, only: t_startf, t_stopf
    use element_mod   , only : Element_t
    real(kind=real_kind)  , intent(in   ) :: edgebuf(nlyr,nbuf)
    integer               , intent(in   ) :: nlyr
    integer               , intent(in   ) :: vlyr
    real(kind=real_kind)  , intent(inout) :: v(np,np,vlyr,tdim,nelemd)
    integer               , intent(in   ) :: kptr
    type(element_t)        ,intent(in   ) :: elem(:)
    integer                ,intent(in   ) :: nets,nete,tdim,tl
    ! Local
    integer :: i,k,ll,is,ie,in,iw,el,kc,kk,glob_k,loc_ind,ii,jj, j
    integer, parameter :: kchunk = 32
    real(kind=real_kind) :: vtmp(np,np,kchunk)
    call t_startf('edge_unpack_max')
    !$acc parallel loop gang collapse(2) present(v,getmapP,edgebuf) private(vtmp)
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
              vtmp(i ,1 ,kk) = max( vtmp(i ,1 ,kk) , edgebuf(kptr+k,getmapP(south,el)+i) )
              vtmp(np,i ,kk) = max( vtmp(np,i ,kk) , edgebuf(kptr+k,getmapP(east ,el)+i) )
              vtmp(i ,np,kk) = max( vtmp(i ,np,kk) , edgebuf(kptr+k,getmapP(north,el)+i) )
              vtmp(1 ,i ,kk) = max( vtmp(1 ,i ,kk) , edgebuf(kptr+k,getmapP(west ,el)+i) )
            endif
          enddo
        enddo
        !$acc loop vector collapse(2)
        do kk = 1 , kchunk
          do i = 1 , max_corner_elem
            k = (kc-1)*kchunk+kk
            if (k <= vlyr) then
              ll = swest+0*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,1 ,kk) = max( vtmp(1 ,1 ,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+1*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,1 ,kk) = max( vtmp(np,1 ,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+2*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(1  ,np,kk) = max( vtmp(1 ,np,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
              ll = swest+3*max_corner_elem+i-1; if(getmapP(ll,el) /= -1) vtmp(np ,np,kk) = max( vtmp(np,np,kk) , edgebuf(kptr+k,getmapP(ll,el)+1) )
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
    call t_stopf('edge_unpack_max')
  end subroutine edgeVunpackMax

end module edge_openacc_mod

