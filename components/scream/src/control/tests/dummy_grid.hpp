#include "share/grid/se_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

namespace scream {

// === A dummy physics grids for this test === //

class DummyPhysicsGrid : public SEGrid
{
public:
  DummyPhysicsGrid (const int num_cols)
   : SEGrid(std::string("Physics"),GridType::SE_NodeBased,num_cols)
  {
    // Nothing to do here
  }

  DummyPhysicsGrid (const int num_cols, const bool fwd)
   : SEGrid(std::string("Physics") + (fwd ? "_fwd" : "_bwd"),GridType::SE_NodeBased,num_cols)
  {
    // Nothing to do here
  }
  ~DummyPhysicsGrid () = default;
};

template<typename ScalarType, typename DeviceType>
class DummyPhysicsGridRemapper : public AbstractRemapper<ScalarType,DeviceType>
{
public:
  using base_type       = AbstractRemapper<ScalarType,DeviceType>;
  using grid_type       = typename base_type::grid_type;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;


  DummyPhysicsGridRemapper(const std::shared_ptr<const grid_type>& grid1,
                           const std::shared_ptr<const grid_type>& grid2)
   : base_type(grid1,grid2)
  {
    const auto& g1n = grid1->name();
    const auto& g2n = grid2->name();

    auto g1 = std::dynamic_pointer_cast<const DummyPhysicsGrid>(grid1);
    auto g2 = std::dynamic_pointer_cast<const DummyPhysicsGrid>(grid2);
    EKAT_REQUIRE_MSG (static_cast<bool>(g1) && static_cast<bool>(g2),
                      "Error! This dummy remapper only works with DummyPhysicsGrid.\n");

    EKAT_REQUIRE_MSG((g1n=="Physics_fwd" && g2n=="Physics_bwd") ||
                     (g1n=="Physics_bwd" && g2n=="Physics_fwd"),
                     "Error! This dummy remapper only works if the two grids are one fwd and the other bwd.\n");
  }

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    // Src and tgt grids are the same, so return the input
    return tgt_layout;
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    // Src and tgt grids are the same, so return the input
    return src_layout;
  }

// Cuda does not allow enclosing function of a lambda to have private or protected access
#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_fields[ifield].first.get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_fields[ifield].second.get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_fields[ifield].first;
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_fields[ifield].second;
  }

  void do_registration_begins () {}
  void do_registration_ends () {}

  void do_register_field (const identifier_type& src, const identifier_type& tgt) override {
    field_type f1(src);
    field_type f2(tgt);
    m_fields.emplace_back(std::make_pair(f1,f2));
  }
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override {
    m_fields[ifield].first  = src;
    m_fields[ifield].second = tgt;
  }
  void do_unregister_field (const int ifield) override {
    m_fields.erase(m_fields.begin()+ifield);
  }

  void do_remap_fwd () const override {
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
  void do_remap_bwd () const override {
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

} // namespace scream
