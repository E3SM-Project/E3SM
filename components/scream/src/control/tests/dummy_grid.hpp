#include "share/grid/point_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

namespace scream {

/*
 * A dummy point grid remapper for testing
 * 
 * The dummy remapper has limited abilities. It can only handle
 * rank 1 and rank 2 fields, from 'Point Grid A' and 'Point Grid B'.
 * Rank 1 fields are simply deep-copied, while rank 2 fields get
 * their layout swapped during remap (meaning a (COL,VL) field
 * is remapped into a (VL,COL) field).
 */

template<typename RealType>
class DummyPointGridRemapper : public AbstractRemapper<RealType>
{
public:
  using base_type       = AbstractRemapper<RealType>;
  using grid_type       = typename base_type::grid_type;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;


  DummyPointGridRemapper(const std::shared_ptr<const grid_type>& grid1,
                         const std::shared_ptr<const grid_type>& grid2)
   : base_type(grid1,grid2)
  {
    const auto& g1n = grid1->name();
    const auto& g2n = grid2->name();

    auto g1 = std::dynamic_pointer_cast<const PointGrid>(grid1);
    auto g2 = std::dynamic_pointer_cast<const PointGrid>(grid2);
    EKAT_REQUIRE_MSG (static_cast<bool>(g1) && static_cast<bool>(g2),
                      "Error! This dummy remapper only works with PointGrid.\n");

    EKAT_REQUIRE_MSG(g1n=="Point Grid A" && g2n=="Point Grid B",
      "Error! This dummy remapper only works if the two grids are called 'Point Grid A' and 'Point Grid B'.\n");
  }

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    const auto rank = tgt_layout.rank();
    EKAT_REQUIRE_MSG (rank==1 || rank==2,
                      "Error! Only rank-1 and rank-2 layouts supported.\n");
    if (rank==1) {
      return tgt_layout;
    }
    const auto& tags = tgt_layout.tags();
    const auto& dims = tgt_layout.dims();
    return FieldLayout( {tags[1],tags[0]}, {dims[1],dims[0]} );
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    const auto rank = src_layout.rank();
    EKAT_REQUIRE_MSG (rank==1 || rank==2,
                      "Error! Only rank-1 and rank-2 layouts supported.\n");
    if (rank==1) {
      return src_layout;
    }
    const auto& tags = src_layout.tags();
    const auto& dims = src_layout.dims();
    return FieldLayout( {tags[1],tags[0]}, {dims[1],dims[0]} );
  }

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const {
    const auto rank = src.rank();
    return rank==tgt.rank() &&
          (rank==1 ? src==tgt :
            (rank==2 && src.tag(0)==tgt.tag(1) && src.tag(1)==tgt.tag(0) && 
                        src.dim(0)==tgt.dim(1) && src.dim(1)==tgt.dim(0)));
  }


protected:

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

// CUDA needs top level lambdas to be enclosed by a method that is public.
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void do_remap_fwd () const override {
    for (const auto& it : m_fields) {
      const auto& layout = it.first.get_header().get_identifier().get_layout();
      if (layout.rank()==1) {
        auto src = it.first.template get_reshaped_view<Real*>();
        auto tgt = it.second.template get_reshaped_view<Real*>();
        Kokkos::deep_copy(tgt,src);
      } else {
        auto src_id = it.first.get_header().get_identifier();
        auto tgt_id = it.second.get_header().get_identifier();

        auto src = it.first.template get_reshaped_view<Real**>();
        auto tgt = it.second.template get_reshaped_view<Real**>();

        auto ncols = this->m_src_grid->get_num_local_dofs();
        auto nlevs = this->m_src_grid->get_num_vertical_levels();
        Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*nlevs),
                             KOKKOS_LAMBDA (int idx) {
          const int icol = idx / nlevs;
          const int ilev = idx % nlevs;
          
          tgt(ilev,icol) = src(icol,ilev);
        });
        Kokkos::fence();
      }
    }
  }
  void do_remap_bwd () const override {
    for (const auto& it : m_fields) {
      const auto& layout = it.first.get_header().get_identifier().get_layout();
      if (layout.rank()==1) {
        auto src = it.first.template get_reshaped_view<Real*>();
        auto tgt = it.second.template get_reshaped_view<Real*>();
        Kokkos::deep_copy(src,tgt);
      } else {
        auto src = it.first.template get_reshaped_view<Real**>();
        auto tgt = it.second.template get_reshaped_view<Real**>();

        auto ncols = this->m_src_grid->get_num_local_dofs();
        auto nlevs = this->m_src_grid->get_num_vertical_levels();
        Kokkos::parallel_for(Kokkos::RangePolicy<>(0,ncols*nlevs),
                             KOKKOS_LAMBDA (int idx) {
          const int icol = idx / nlevs;
          const int ilev = idx % nlevs;
          
          src(icol,ilev) = tgt(ilev,icol);
        });
        Kokkos::fence();
      }
    }
  }

protected:
  std::vector<std::pair<field_type,field_type>> m_fields;
};

} // namespace scream
