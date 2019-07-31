#include "share/grid/default_grid.hpp"
#include "share/remap/abstract_remapper.hpp"

namespace scream {

// === A dummy physics grids for this test === //

class DummyPhysicsGrid : public DefaultGrid<GridType::Physics>
{
public:
  DummyPhysicsGrid (const int num_cols, const bool fwd)
   : DefaultGrid<GridType::Physics>(std::string("Physics") + (fwd ? "_fwd" : "_bwd"))
  {
    m_num_dofs = num_cols;

  }
  ~DummyPhysicsGrid () = default;

protected:
};

template<typename ScalarType, typename DeviceType>
class DummyPhysicsGridRemapper : public AbstractRemapper<ScalarType,DeviceType>
{
public:
  using base_type       = AbstractRemapper<ScalarType,DeviceType>;
  using grid_type       = typename base_type::grid_type;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;


  DummyPhysicsGridRemapper(const std::shared_ptr<const grid_type>& grid1,
                           const std::shared_ptr<const grid_type>& grid2)
   : base_type(grid1,grid2)
  {
    const auto& g1n = grid1->name();
    const auto& g2n = grid2->name();

    auto g1 = std::dynamic_pointer_cast<const DummyPhysicsGrid>(grid1);
    auto g2 = std::dynamic_pointer_cast<const DummyPhysicsGrid>(grid2);
    scream_require_msg (static_cast<bool>(g1) && static_cast<bool>(g2),
                        "Error! This dummy remapper only works with DummyPhysicsGrid.\n");

    scream_require_msg((g1n=="Physics_fwd" && g2n=="Physics_bwd") ||
                       (g1n=="Physics_bwd" && g2n=="Physics_fwd"),
                       "Error! This dummy remapper only works if the two grids are one fwd and the other bwd.\n");
  }

  base_type* clone () const {
    auto remapper = new DummyPhysicsGridRemapper(this->get_src_grid(),this->get_tgt_grid());

    remapper->m_fields.reserve(this->m_fields.size());
    for (auto& it : this->m_fields) {
      remapper->m_fields.emplace_back(it);
    }

    return remapper;
  }

  // One grid is backward of the other, but same layouts
  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const {
    return tgt_layout;
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const {
    return src_layout;
  }

protected:
  const identifier_type& do_get_src_field_id (const int ifield) const {
    return m_fields[ifield].first.get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const {
    return m_fields[ifield].second.get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const {
    return m_fields[ifield].first;
  }
  const field_type& do_get_tgt_field (const int ifield) const {
    return m_fields[ifield].second;
  }

  void do_registration_begins () {}
  void do_registration_complete () {}

  void do_register_field (const identifier_type& src, const identifier_type& tgt) {
    field_type f1(src);
    field_type f2(tgt);
    m_fields.emplace_back(std::make_pair(f1,f2));
  }
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) {
    m_fields[ifield].first  = src;
    m_fields[ifield].second = tgt;
  }

  void do_remap_fwd () const {
    for (int i=0; i<this->m_num_fields; ++i) {
      auto src = m_fields[i].first.get_view();
      auto tgt = m_fields[i].second.get_view();

      auto n = src.size();
      Kokkos::parallel_for(Kokkos::RangePolicy<>(0,n),
                           KOKKOS_LAMBDA (int i) {
        tgt[i] = src[n-i-1];
      });
      Kokkos::fence();
    }
  }
  void do_remap_bwd () const {
    for (int i=0; i<this->m_num_fields; ++i) {
      auto src = m_fields[i].first.get_view();
      auto tgt = m_fields[i].second.get_view();

      auto n = src.size();
      Kokkos::parallel_for(Kokkos::RangePolicy<>(0,n),
                           KOKKOS_LAMBDA (int i) {
        src[i] = tgt[n-i-1];
      });
      Kokkos::fence();
    }
  }

  std::vector<std::pair<field_type,field_type>> m_fields;
};

} // empty namespace
