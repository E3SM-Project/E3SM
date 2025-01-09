#ifndef INCLUDE_COMPOSE_CEDR_CDR_HPP
#define INCLUDE_COMPOSE_CEDR_CDR_HPP

#include "compose_cedr.hpp"
#include "compose_cedr_qlt.hpp"
#include "compose_cedr_caas.hpp"

namespace homme {

struct Alg {
  enum Enum { qlt, qlt_super_level, qlt_super_level_local_caas, caas,
              caas_super_level };
  static Enum convert (int cdr_alg) {
    switch (cdr_alg) {
    case 2:  return qlt;
    case 20: return qlt_super_level;
    case 21: return qlt_super_level_local_caas;
    case 3:  return caas;
    case 30: return caas_super_level;
    case 42: return caas_super_level; // actually none
    default: cedr_throw_if(true,  "cdr_alg " << cdr_alg << " is invalid.");
    }
  }
  static bool is_qlt (Enum e) {
    return (e == qlt || e == qlt_super_level ||
            e == qlt_super_level_local_caas);
  }
  static bool is_caas (Enum e) {
    return e == caas || e == caas_super_level;
  }
  static bool is_point (Enum e) {
    return false;
  }
  static bool is_suplev (Enum e) {
    return (e == qlt_super_level || e == caas_super_level ||
            e == qlt_super_level_local_caas);
  }
};

template <typename MT>
struct CDR {
  typedef cedr::Int Int;
  typedef cedr::Real Real;

  typedef std::shared_ptr<CDR> Ptr;
  typedef compose::QLT<typename MT::DES> QLTT;
#ifdef COMPOSE_PORT
  typedef cedr::caas::CAAS<typename MT::DES> CAAST;
#else
  typedef compose::CAAS CAAST;
#endif

  typedef Kokkos::View<Int*, typename MT::DES> Idxs;
  typedef typename Idxs::HostMirror IdxsH;

  typedef Kokkos::View<char*, typename MT::DES> Bools;
  typedef typename Bools::HostMirror BoolsH;

  enum { nsublev_per_suplev = 8 };
  
  const Alg::Enum alg;
  const Int ncell, nlclcell, nlev, np, qsize, nsublev, nsuplev;
  const bool threed, cdr_over_super_levels, caas_in_suplev, hard_zero;
  const cedr::mpi::Parallel::Ptr p;
  cedr::tree::Node::Ptr tree; // Don't need this except for unit testing.
  cedr::CDR::Ptr cdr;
  Idxs ie2gci; // Map Homme ie to Homme global cell index.
  Idxs ie2lci; // Map Homme ie to CDR local cell index (lclcellidx).
  IdxsH ie2lci_h, ie2gci_h;
  Bools nonneg;
  BoolsH nonneg_h;
  bool run; // for debugging, it can be useful not to run the CEDR.

  CDR(Int cdr_alg_, Int ngblcell_, Int nlclcell_, Int nlev_, Int np_, Int qsize_,
      bool use_sgi, bool independent_time_steps, const bool hard_zero_,
      const Int* gid_data, const Int* rank_data, const cedr::mpi::Parallel::Ptr& p_,
      Int fcomm);

  CDR(const CDR&) = delete;
  CDR& operator=(const CDR&) = delete;

  void init_tracers(const bool need_conservation);

  void get_buffers_sizes(size_t& s1, size_t &s2);

  void set_buffers(Real* b1, Real* b2);

private:
  bool inited_tracers_;
};

} // namespace homme

#endif
