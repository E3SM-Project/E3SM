#ifndef INCLUDE_COMPOSE_SLMM_ISLMPI_HPP
#define INCLUDE_COMPOSE_SLMM_ISLMPI_HPP

#include "compose.hpp"
#include "compose_homme.hpp"
#include "compose_slmm.hpp"
#include "compose_slmm_advecter.hpp"

#include <mpi.h>

#include <memory>

// AMB 2017/06-2020/05 Initial for E3SMv2
// AMB 2020/05-2021/01 Performance-portable impl
// AMB 2021/04         Support doubly-periodic planar mode
// AMB 2024/04-2025/01 Enhanced trajectory method

namespace homme {
namespace mpi { //todo Share with cedr.

class Parallel {
  MPI_Comm comm_;
public:
  typedef std::shared_ptr<Parallel> Ptr;
  Parallel (MPI_Comm comm) : comm_(comm) {}
  MPI_Comm comm () const { return comm_; }
  Int size () const {
    int sz = 0;
    MPI_Comm_size(comm_, &sz);
    return sz;
  }
  Int rank () const {
    int pid = 0;
    MPI_Comm_rank(comm_, &pid);
    return pid;
  }
  Int root () const { return 0; }
  bool amroot () const { return rank() == root(); }
};

inline Parallel::Ptr make_parallel (MPI_Comm comm) {
  return std::make_shared<Parallel>(comm);
}

struct Request {
  MPI_Request request;

#ifdef COMPOSE_DEBUG_MPI
  int unfreed;
  Request();
  ~Request();
#endif
};

template <typename T> MPI_Datatype get_type();
template <> inline MPI_Datatype get_type<int>() { return MPI_INT; }
template <> inline MPI_Datatype get_type<double>() { return MPI_DOUBLE; }
template <> inline MPI_Datatype get_type<long>() { return MPI_LONG_INT; }

template <typename T>
int isend (const Parallel& p, const T* buf, int count, int dest, int tag,
           Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Isend(const_cast<T*>(buf), count, dt, dest, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
}

template <typename T>
int irecv (const Parallel& p, T* buf, int count, int src, int tag, Request* ireq) {
  MPI_Datatype dt = get_type<T>();
  MPI_Request ureq;
  MPI_Request* req = ireq ? &ireq->request : &ureq;
  int ret = MPI_Irecv(buf, count, dt, src, tag, p.comm(), req);
  if ( ! ireq) MPI_Request_free(req);
#ifdef COMPOSE_DEBUG_MPI
  else ireq->unfreed++;
#endif
  return ret;
}

int waitany(int count, Request* reqs, int* index, MPI_Status* stats = nullptr);
int waitall(int count, Request* reqs, MPI_Status* stats = nullptr);
int wait(Request* req, MPI_Status* stat = nullptr);

template <typename T>
int all_reduce (const Parallel& p, const T* sendbuf, T* rcvbuf, int count, MPI_Op op) {
  MPI_Datatype dt = get_type<T>();
  return MPI_Allreduce(const_cast<T*>(sendbuf), rcvbuf, count, dt, op, p.comm());
}
} // namespace mpi

namespace islmpi {
// Meta and bulk data for the interpolation SL method with special comm pattern.

#ifdef COMPOSE_BOUNDS_CHECK
# define slmm_assert_high(condition) slmm_assert(condition)
# define slmm_kernel_assert_high(condition) slmm_kernel_assert(condition)
#else
# define slmm_assert_high(condition)
# define slmm_kernel_assert_high(condition)
#endif

// FixedCapList, ListOfLists, and BufferLayoutArray are simple and somewhat
// problem-specific array data structures for use in IslMpi.
template <typename T, typename DT>
struct FixedCapList {
  typedef ko::View<T*, DT> Array;
  typedef FixedCapList<T, typename Array::host_mirror_space> Mirror;

  // Empty ctor b/c we need to permit view of view construction.
  SLMM_KIF FixedCapList () {}

  FixedCapList (const Int& cap) {
    slmm_assert_high(cap >= 0);
    reset_capacity(cap);
  }

  void reset_capacity (const Int& cap, const bool also_size = false) {
    slmm_assert(cap >= 0);
    if (d_.size() == 0) init_n();
    ko::resize(d_, cap);
    set_n_from_host(also_size ? cap : 0);
  }

  SLMM_KIF Int capacity () const { return d_.size(); }

  // If empty ctor was called, nothing below here is valid until reset_capacity
  // is called.

  SLMM_KIF void clear () const { set_n(0); }
  SLMM_KIF Int size () const { return n(); }
  SLMM_KIF T& operator() (const Int& i) const { slmm_kernel_assert_high(i >= 0 && i < n()); return d_[i]; }
  SLMM_KIF void inc () const { ++get_n_ref(); slmm_kernel_assert_high(n() <= static_cast<Int>(d_.size())); }
  SLMM_KIF void inc (const Int& dn) const { get_n_ref() += dn; slmm_kernel_assert_high(n() <= static_cast<Int>(d_.size())); }
  SLMM_KIF T& atomic_inc_and_return_next () const {
#ifdef COMPOSE_PORT
    volatile Int* n_vol = &n_();
    Int n_read = -1;
    for (;;) {
      n_read = *n_vol;
      if (ko::atomic_compare_exchange_strong(n_vol, n_read, n_read+1))
        break;
    }
    return d_[n_read];
#else
    inc();
    return back();
#endif
  }

  SLMM_KIF T* data () const { return d_.data(); }  
  SLMM_KIF T& back () const { slmm_kernel_assert_high(n() > 0); return d_[n()-1]; }
  SLMM_KIF T* begin () const { return d_.data(); }
  SLMM_KIF T* end () const { return d_.data() + n(); }

  void zero () { if (d_.size()) ko::deep_copy(d_, 0); }

  // Copy from s to this.
  template <typename DTS>
  void copy (const FixedCapList<T, DTS>& s) {
    siqk::resize_and_copy(d_, s.view());
    init_n();
    set_n(s);
  }

  // Use everything that follows only for low-level host-device things.

  SLMM_KIF void set_n (const Int& n0) const { get_n_ref() = n0; }

  template <typename DTS>
  void set_n (const DTS& s) const {
#ifdef COMPOSE_PORT
    ko::deep_copy(n_, s.n());
#else
    get_n_ref() = s.n();
#endif
  }

#ifdef COMPOSE_PORT
  typedef ko::View<Int, DT> NT;

  SLMM_KIF Int n () const { return n_(); }

  void init_n () { n_ = NT("FixedCapList::n_"); }
  void set_n_from_host (const Int& n0) { ko::deep_copy(n_, n0); }

  // Create a FixedCapList whose View is a mirror view of this.
  Mirror mirror () const {
    Mirror v;
    v.set_view(ko::create_mirror_view(d_));
    v.set_n_view(ko::create_mirror_view(n_));
    ko::deep_copy(v.n_view(), n_);
    return v;
  }

  const ko::View<Int, DT>& n_view () const { return n_; }
  void set_n_view (const NT& v) { n_ = v; }
#else
  typedef Int NT;

  SLMM_KIF Int n () const { return n_; }

  void init_n () {}
  void set_n_from_host (const Int& n0) { set_n(n0); }

  Mirror mirror () const {
    Mirror v;
    v.set_view(ko::create_mirror_view(d_));
    v.set_n(n_);
    return v;
  }
#endif

  const Array& view () const { return d_; }
  void set_view (const Array v) { d_ = v; }

  FixedCapList<T,DT> unmanaged () const {
    FixedCapList<T,DT> ufcl;
    ufcl.d_ = Array(d_.data(), d_.extent(0));
#ifdef COMPOSE_PORT
    ufcl.n_ = ko::View<Int, DT>(n_.data());
#else
    ufcl.n_ = n_;
#endif
    return ufcl;
  }

private:
  Array d_;

#ifndef COMPOSE_PORT
  // You'll notice in a number of spots that there is strange const/mutable
  // stuff going on. This is driven by what Kokkos needs.
  mutable
#endif
  NT n_;

#ifdef COMPOSE_PORT
  SLMM_KIF Int& get_n_ref () const { return n_(); }
#else
  SLMM_KIF Int& get_n_ref () const { return n_; }
#endif
};

template <typename T, typename DTD, typename DTS>
void deep_copy (FixedCapList<T, DTD>& d, const FixedCapList<T, DTS>& s) {
  slmm_assert_high(d.capacity() == s.capacity());
  if (d.view().size() > 0) ko::deep_copy(d.view(), s.view());
#ifdef COMPOSE_PORT
  ko::deep_copy(d.n_view(), s.n_view());
#else
  d.set_n(s.n());
#endif
}

template <typename T>
struct FixedCapListHostOnly {
  FixedCapListHostOnly (const Int cap = 0) {
    slmm_assert_high(cap >= 0);
    reset_capacity(cap);
  }

  void reset_capacity (const Int cap, const bool also_size = false) {
    slmm_assert(cap >= 0);
    d_.resize(cap);
    n_ = also_size ? cap : 0;
  }

  Int capacity () const { return d_.size(); }
  Int size () const { return n_; }
  Int n () const { return n_; }
  
  void clear () { n_ = 0; }

  void inc () { ++n_; slmm_kernel_assert_high(n_ <= static_cast<Int>(d_.size())); }
  void inc (const Int& dn) { n_ += dn; slmm_kernel_assert_high(n_ <= static_cast<Int>(d_.size())); }

  T& operator() (const Int& i) { slmm_kernel_assert_high(i >= 0 && i < n_); return d_[i]; }

  T* data () { return d_.data(); }  
  T& back () { slmm_kernel_assert_high(n_ > 0); return d_[n_-1]; }
  T* begin () { return d_.data(); }
  T* end () { return d_.data() + n_; }

private:
  std::vector<T> d_;
  Int n_;
};

template <typename DT> struct BufferLayoutArray;

template <typename T, typename DT>
struct ListOfLists {
  template <typename T1> using Array = ko::View<T1*, DT>;
  typedef ListOfLists<T, typename Array<T>::host_mirror_space> Mirror;

  struct List {
    SLMM_KIF Int n () const { return n_; }

    SLMM_KIF T& operator() (const Int& i) const {
      slmm_kernel_assert_high(i >= 0 && i < n_); return d_[i];
    }

    SLMM_KIF T* data () const { return d_; }
    SLMM_KIF T* begin () const { return d_; }
    SLMM_KIF T* end () const { return d_ + n_; }
    SLMM_KIF Array<T> view () const { return Array<T>(d_, n_); }

  private:
    friend class ListOfLists<T, DT>;
    SLMM_KIF List (T* d, const Int& n) : d_(d), n_(n) { slmm_kernel_assert_high(n_ >= 0); }
    T* const d_;
    const Int n_;
  };

  ListOfLists () {}
  ListOfLists (const Int nlist, const Int* nlist_per_list) { init(nlist, nlist_per_list); }
  void init (const Int nlist, const Int* nlist_per_list, T* buf = nullptr) {
    slmm_assert(nlist >= 0);
    ptr_ = Array<Int>("ptr_", nlist+1);
    ptr_h_ = ko::create_mirror_view(ptr_);
    ptr_h_[0] = 0;
    for (Int i = 0; i < nlist; ++i) {
      slmm_assert(nlist_per_list[i] >= 0);
      ptr_h_[i+1] = ptr_h_[i] + nlist_per_list[i];
    }
    ko::deep_copy(ptr_, ptr_h_);
    if (buf) {
      d_ = Array<T>(buf, ptr_h_[nlist]);
    } else {
      d_ = Array<T>("d_", ptr_h_[nlist]);
    }
  }

  SLMM_KIF Int n () const { return static_cast<Int>(ptr_.size()) - 1; }
  SLMM_KIF List operator() (const Int& i) const {
    slmm_kernel_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1);
    return List(const_cast<T*>(&d_[ptr_[i]]), ptr_[i+1] - ptr_[i]);
  }
  List get_h (const Int& i) const {
    slmm_kernel_assert_high(i >= 0 && i < static_cast<Int>(ptr_h_.size()) - 1);
    return List(const_cast<T*>(d_.data() + ptr_h_[i]), ptr_h_[i+1] - ptr_h_[i]);
  }
  SLMM_KIF T& operator() (const Int& i, const Int& j) const {
    slmm_kernel_assert_high(i >= 0 && i < static_cast<Int>(ptr_.size()) - 1 &&
                            j >= 0 && j < ptr_[i+1] - ptr_[i]);
    return d_[ptr_[i] + j];
  }

  void zero () {
#ifdef COMPOSE_HORIZ_OPENMP
#   pragma omp for
    for (Int i = 0; i < ptr_[n()]; ++i)
      d_[i] = 0;
#else
    ko::deep_copy(d_, 0);
#endif
  }

  SLMM_KIF T* data () const { return d_.data(); }

  ListOfLists<T,DT> unmanaged () const {
    ListOfLists<T,DT> ulol;
    ulol.d_ = Array<T>(d_.data(), d_.extent(0));
    ulol.ptr_ = Array<Int>(ptr_.data(), ptr_.extent(0));
    ulol.ptr_h_ = typename Array<Int>::HostMirror(ptr_h_.data(), ptr_h_.extent(0));
    return ulol;
  }

  // For device-host stuff:

  void set_views (const Array<T>& d, const Array<Int>& ptr,
                  const typename Array<Int>::HostMirror& ptr_h) {
    d_ = d; ptr_ = ptr; ptr_h_ = ptr_h;
  }

  Mirror mirror () const {
    Mirror v;
    const auto ptr = ko::create_mirror_view(ptr_);
    ko::deep_copy(ptr, ptr_h_);
    v.set_views(ko::create_mirror_view(d_), ptr, ptr_h_);
    return v;
  }

  const Array<T>& d_view () const { return d_; }
  const Array<Int>& ptr_view () const { return ptr_; }
  const typename Array<Int>::HostMirror& ptr_h_view () const { return ptr_h_; }

private:
  friend class BufferLayoutArray<DT>;
  Array<T> d_;
  Array<Int> ptr_;
  typename Array<Int>::HostMirror ptr_h_;
};

struct LayoutTriple {
  Int xptr, qptr, cnt;
  SLMM_KIF LayoutTriple () : LayoutTriple(0) {}
  SLMM_KIF LayoutTriple (const Int& val) { xptr = qptr = cnt = 0; }
};

template <typename T, typename DTD, typename DTS>
void deep_copy (ListOfLists<T, DTD>& d, const ListOfLists<T, DTS>& s) {
  ko::deep_copy(d.d_view(), s.d_view());
  ko::deep_copy(d.ptr_view(), s.ptr_view());
  ko::deep_copy(d.ptr_h_view(), s.ptr_h_view());
}

template <typename DT>
struct BufferLayoutArray {
  struct BufferRankLayoutArray {
    SLMM_KIF LayoutTriple& operator() (const Int& lidi, const Int& lev) const {
      slmm_kernel_assert_high(lidi >= 0 && lev >= 0 && lidi*nlev_ + lev < d_.n());
      return d_(lidi*nlev_ + lev);
    }

  private:
    friend class BufferLayoutArray;
    SLMM_KIF BufferRankLayoutArray (const typename ListOfLists<LayoutTriple, DT>::List& d,
                                    const Int& nlev)
      : d_(d), nlev_(nlev) {}
    typename ListOfLists<LayoutTriple, DT>::List d_;
    Int nlev_;
  };

  BufferLayoutArray () : nlev_(0) {}
  BufferLayoutArray (const Int& nrank, const Int* nlid_per_rank, const Int& nlev)
  { init(nrank, nlid_per_rank, nlev); }

  void init (const Int& nrank, const Int* nlid_per_rank, const Int& nlev) {
    slmm_assert(nrank >= 0 && nlev > 0);
    nlev_ = nlev;
    std::vector<Int> ns(nrank);
    for (Int i = 0; i < nrank; ++i) {
      slmm_assert(nlid_per_rank[i] > 0);
      ns[i] = nlid_per_rank[i] * nlev;
    }
    d_.init(nrank, ns.data());
  }

  Int nlev () const { return nlev_; }

  void zero () { d_.zero(); }

  SLMM_KIF LayoutTriple& operator() (const Int& ri, const Int& lidi, const Int& lev) const {
    slmm_kernel_assert_high(ri >= 0 && ri < d_.n() &&
                            lidi >= 0 && lev >= 0 &&
                            lidi*nlev_ + lev < d_(ri).n());
    return d_.data()[d_.ptr_[ri] + lidi*nlev_ + lev];
  }

  SLMM_KIF BufferRankLayoutArray operator() (const Int& ri) const {
    slmm_kernel_assert_high(ri >= 0 && ri < d_.n());
    return BufferRankLayoutArray(d_(ri), nlev_);
  }

  // For device-host stuff:

  typedef BufferLayoutArray<typename ko::View<Int,DT>::host_mirror_space> Mirror;

  ListOfLists<LayoutTriple, DT>& get_lol () { return d_; }
  const ListOfLists<LayoutTriple, DT>& get_lol () const { return d_; }
  void set_lol (const ListOfLists<LayoutTriple, DT>& d) { d_ = d; }
  void set_nlev (const Int& nlev) { nlev_ = nlev; }

  Mirror mirror () const {
    Mirror v;
    const auto m = d_.mirror();
    v.set_lol(m);
    v.set_nlev(nlev_);
    return v;
  }

  BufferLayoutArray<DT> unmanaged () const {
    BufferLayoutArray<DT> ubla;
    ubla.d_ = d_.unmanaged();
    ubla.nlev_ = nlev_;
    return ubla;
  }

private:
  ListOfLists<LayoutTriple, DT> d_;
  Int nlev_;
};

template <typename DTD, typename DTS>
void deep_copy (BufferLayoutArray<DTD>& d, const BufferLayoutArray<DTS>& s) {
  slmm_assert(d.nlev() == s.nlev());
  deep_copy(d.get_lol(), s.get_lol());
}

struct GidRank {
  Int
    gid,      // cell global ID
    rank,     // the rank that owns the cell
    rank_idx, // index into list of ranks with whom I communicate, including me
    lid_on_rank,     // the local ID of the cell on the owning rank
    lid_on_rank_idx; // index into list of LIDs on the rank
};
struct OwnItem {
  short lev;   // level index
  short k;     // linearized GLL index
};
struct RemoteItem {
  Int q_extrema_ptr, q_ptr; // pointers into recvbuf
  short lev, k;
};

// Meta and bulk data for the interpolation SL MPI communication pattern.
template <typename MT = ko::MachineTraits>
struct IslMpi {
  typedef typename MT::HES HES;
  typedef typename MT::DES DES;
  typedef typename MT::HDT HDT;
  typedef typename MT::DDT DDT;

  using Advecter = slmm::Advecter<MT>;

  typedef std::shared_ptr<IslMpi> Ptr;

  template <typename Datatype, typename DT>
  using Array = ko::View<Datatype, siqk::Layout, DT>;
  template <typename Datatype>
  using ArrayH = ko::View<Datatype, siqk::Layout, HDT>;
  template <typename Datatype>
  using ArrayD = ko::View<Datatype, siqk::Layout, DDT>;

  // The comm and real data associated with an element patch, the set of
  // elements surrounding an owned cell.
  template <typename DT>
  struct ElemData {
    GidRank* me;                      // the owned cell
    FixedCapList<GidRank, DT> nbrs;   // the cell's neighbors (but including me)
    Int nin1halo;                     // nbrs[0:n]
    FixedCapList<OwnItem, DT> own;    // points whose q are computed with own rank's data
    FixedCapList<RemoteItem, DT> rmt; // points computed by a remote rank's data
    Array<Int**, DT> src;             // src(lev,k) = get_src_cell
    Array<Real**[2], DT> q_extrema;
    const Real* qdp, * dp;  // the owned cell's data
    Real* q;
  };

  typedef ElemData<HDT> ElemDataH;
  typedef ElemData<DDT> ElemDataD;
  typedef FixedCapList<ElemDataH, HDT> ElemDataListH;
  typedef FixedCapList<ElemDataD, DDT> ElemDataListD;

  const mpi::Parallel::Ptr p;
  const typename Advecter::ConstPtr advecter;
  const Int np, np2, nlev, qsize, qsized, nelemd, halo;
  const bool traj_3d;
  const Int traj_nsubstep, dep_points_ndim;

  Real etai_beg, etai_end;
  ArrayD<Real*> etam;

  ElemDataListH ed_h; // this rank's owned cells, indexed by LID
  ElemDataListD ed_d;
  typename ElemDataListD::Mirror ed_m; // handle managed allocs

  const typename TracerArrays<MT>::Ptr tracer_arrays;

  // IDs.
  FixedCapList<Int, HDT> ranks, mylid_with_comm_tid_ptr_h;
  FixedCapList<Int, DDT> nx_in_rank, mylid_with_comm_d;
  ListOfLists <Int, DDT> nx_in_lid, lid_on_rank;
  BufferLayoutArray<DDT> bla;

  // MPI comm data.
  FixedCapListHostOnly<mpi::Request> sendreq, recvreq;
  FixedCapList<Int, HDT> recvreq_ri;
  ListOfLists<Real, DDT> sendbuf, recvbuf;
#ifdef COMPOSE_MPI_ON_HOST
  typename ListOfLists<Real, DDT>::Mirror sendbuf_h, recvbuf_h;
#endif
  FixedCapList<Int, DDT> sendcount, x_bulkdata_offset;
  ListOfLists<Real, HDT> sendbuf_meta_h, recvbuf_meta_h; // not mirrors
  FixedCapList<Int, DDT> rmt_xs, rmt_qs_extrema;
  Int nrmt_xs, nrmt_qs_extrema;

  // Mirror views.
  typename FixedCapList<Int, DDT>::Mirror nx_in_rank_h, sendcount_h,
    x_bulkdata_offset_h, rmt_xs_h, rmt_qs_extrema_h, mylid_with_comm_h;
  typename ListOfLists <Int, DDT>::Mirror nx_in_lid_h, lid_on_rank_h;
  typename BufferLayoutArray<DDT>::Mirror bla_h;

  bool horiz_openmp;
#ifdef COMPOSE_HORIZ_OPENMP
  ListOfLists<omp_lock_t, HDT> ri_lidi_locks;
#endif

  // temporary work space
  std::vector<Int> nlid_per_rank, sendsz, recvsz, sendmetasz, recvmetasz;
  ArrayD<Real**> rwork;

  typedef ArrayD<char***> DepMask;
  typedef ArrayD<Int*[3]> DepList;
  DepMask own_dep_mask;
  DepList own_dep_list;
  Int own_dep_list_len;

  IslMpi (const mpi::Parallel::Ptr& ip, const typename Advecter::ConstPtr& advecter,
          const typename TracerArrays<MT>::Ptr& itracer_arrays,
          Int inp, Int inlev, Int iqsize, Int iqsized, Int inelemd, Int ihalo,
          Int itraj_3d, Int itraj_nsubstep)
    : p(ip), advecter(advecter),
      np(inp), np2(np*np), nlev(inlev), qsize(iqsize), qsized(iqsized), nelemd(inelemd),
      halo(ihalo), traj_3d(itraj_3d), traj_nsubstep(itraj_nsubstep),
      dep_points_ndim(traj_3d && traj_nsubstep > 0 ? 4 : 3),
      tracer_arrays(itracer_arrays)
  {}

  IslMpi(const IslMpi&) = delete;
  IslMpi& operator=(const IslMpi&) = delete;

  ~IslMpi () {
#ifdef COMPOSE_HORIZ_OPENMP
    const Int nrmtrank = static_cast<Int>(ranks.n()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      auto&& locks = ri_lidi_locks(ri);
      for (auto& lock: locks) {
        // The following call is causing a seg fault on at least one Cray KNL
        // system. It doesn't kill the run b/c it occurs after main exits. Not
        // calling this is a memory leak, but it's innocuous in this case. Thus,
        // comment it out:
        //omp_destroy_lock(&lock);
        // It's possible that something in the OpenMP runtime shuts down before
        // this call, and that's the problem. If that turns out to be it, I
        // should make a compose_finalize() call somewhere.
      }
    }
#endif
  }
};

inline int get_tid () {
#ifdef COMPOSE_HORIZ_OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

inline int get_num_threads () {
#ifdef COMPOSE_HORIZ_OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

template <typename MT>
void alloc_mpi_buffers (IslMpi<MT>& cm, Real* sendbuf = nullptr, Real* recvbuf = nullptr);

template <typename MT>
void setup_comm_pattern(IslMpi<MT>& cm, const Int* nbr_id_rank, const Int* nirptr);

namespace extend_halo {
template <typename MT>
void extend_local_meshes(const mpi::Parallel& p,
                         const typename IslMpi<MT>::ElemDataListH& eds,
                         typename IslMpi<MT>::Advecter& advecter);
} // namespace extend_halo

template <typename MT>
void analyze_dep_points(IslMpi<MT>& cm, const Int& nets, const Int& nete,
                        const DepPoints<MT>& dep_points);

template <typename MT>
void init_mylid_with_comm_threaded(IslMpi<MT>& cm, const Int& nets, const Int& nete);
template <typename MT>
void setup_irecv(IslMpi<MT>& cm, const bool skip_if_empty = false);
template <typename MT>
void isend(IslMpi<MT>& cm, const bool want_req = true, const bool skip_if_empty = false);
template <typename MT>
void recv_and_wait_on_send(IslMpi<MT>& cm);
template <typename MT>
void wait_on_send (IslMpi<MT>& cm, const bool skip_if_empty = false);
template <typename MT>
void recv(IslMpi<MT>& cm, const bool skip_if_empty = false);

template <typename MT>
void pack_dep_points_sendbuf_pass1(IslMpi<MT>& cm, const bool trajectory = false);
template <typename MT>
void pack_dep_points_sendbuf_pass2(IslMpi<MT>& cm, const DepPoints<MT>& dep_points,
                                   const bool trajectory = false);

template <typename MT>
void calc_q_extrema(IslMpi<MT>& cm, const Int& nets, const Int& nete);

template <typename MT>
void calc_rmt_q_pass1(IslMpi<MT>& cm, const bool trajectory = false);
template <typename MT>
void calc_rmt_q(IslMpi<MT>& cm);
template <typename MT>
void calc_own_q(IslMpi<MT>& cm, const Int& nets, const Int& nete,
                const DepPoints<MT>& dep_points,
                const QExtrema<MT>& q_min, const QExtrema<MT>& q_max);
template <typename MT>
void copy_q(IslMpi<MT>& cm, const Int& nets,
            const QExtrema<MT>& q_min, const QExtrema<MT>& q_max);

/* Take a semi-Lagrangian step, excluding property preservation.
     dep_points is const in principle, but if
       lev <= semi_lagrange_nearest_point_lev,
   a departure point may be altered if the winds take it outside of the comm
   halo.
     This code was originally developed for the unit sphere. Later, we modified
   it to handle a doubly-periodic plane || to the x-y plane. As a result, in a
   number of spots in the code, one sees the word 'sphere', such as
   SphereToRef::calc_sphere_to_ref, but in fact the class and routine may handle
   the plane, as well.
     For the plane, one must account for periodicity in the coordinate values.
   We refer to 'periodic' and 'continuous' values. Periodic values are those
   that are within the plane's bounds. Continuous can be out of bounds but,
   together with other coordinate values in a patch, form a continuous function.
   LocalMesh::p is initialized to have continuous coordinate values anchored at
   LocalMesh::tgt_elem, which itself has periodic values.
     The departure points that are input to 'step' should have continuous
   coordinate values. That is, if a departure point advects outside of the
   plane's bounds, the values should stay as such. Impl reason: get_src_cell
   operates using LocalMesh, so it naturally wants the departure points to have
   continuous values to match the patch values.
     Once that calculation is done, dep_points is changed to periodic values.
*/
template <typename MT = ko::MachineTraits>
void step(
  IslMpi<MT>& cm, const Int nets, const Int nete,
  Real* dep_points_r,
  Real* q_min_r, Real* q_max_r);

template <typename MT = ko::MachineTraits>
void set_hvcoord(IslMpi<MT>& cm, const Real etai_beg, const Real etai_end,
                 const Real* etam);

template <typename MT = ko::MachineTraits>
void calc_v_departure(
  IslMpi<MT>& cm, const Int nets, const Int nete, const Int step, const Real dtsub,
  Real* dep_points_r, const Real* vnode, Real* vdep);

} // namespace islmpi
} // namespace homme

#endif
