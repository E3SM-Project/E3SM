#ifndef INCLUDE_COMPOSE_CEDR_QLT_HPP
#define INCLUDE_COMPOSE_CEDR_QLT_HPP

#include "cedr_qlt.hpp"

namespace homme {
namespace compose {

class VerticalLevelsData;

template <typename ES>
class QLT : public cedr::qlt::QLT<ES> {
  typedef cedr::qlt::QLT<ES> Super;
  typedef cedr::Int Int;

  std::shared_ptr<VerticalLevelsData> vld_;

  void reconcile_vertical(const Int problem_type, const Int bd_os,
                          const Int bis, const Int bie);
  void runimpl();

public:
  QLT(const cedr::mpi::Parallel::Ptr& p, const cedr::Int& ncells,
      const cedr::tree::Node::Ptr& tree, const cedr::CDR::Options& options,
      const cedr::Int& vertical_levels);

  void run() override;

  static Int unittest();
};

} // namespace compose
} // namespace homme

#endif
