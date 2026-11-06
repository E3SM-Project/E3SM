#include "vertical_level_index.hpp"
#include "level_index_utils.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_units.hpp>

namespace
{
// Find the first position in the array x[0,n) that is past v, that is, with x>v if the
// array is increasing, or with x<v if it is decreasing. If there is none, return n.
// The array must be monotonic (in the direction specified).
KOKKOS_INLINE_FUNCTION
int first_past (const scream::Real* x, const int n, const scream::Real v, const bool descending)
{
  int lo = 0;
  int hi = n;
  while (lo<hi) {
    const int mid = lo + (hi-lo)/2;
    const bool past = descending ? (x[mid]<v) : (x[mid]>v);
    if (past) {
      hi = mid;
    } else {
      lo = mid+1;
    }
  }
  return lo;
}
} // anonymous namespace

namespace scream
{

VerticalLevelIndex::
VerticalLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  const level_index_utils::VerticalTarget target (m_params.get<std::string>("vertical_coordinate"),
                                                  m_params.get<std::string>("value"),
                                                  m_params.get<std::string>("units"));
  m_target = target.target;
  m_descending = target.descending;
  m_extrapolate_bottom = target.extrapolate_bottom;

  const auto& layer = m_params.get<std::string>("vertical_layer");
  EKAT_REQUIRE_MSG (layer=="mid" or layer=="int",
      "Error! Invalid vertical layer for VerticalLevelIndex.\n"
      " - input value : " + layer + "\n"
      " - valid values: 'mid', 'int'\n");

  m_coord_name = target.coord_field(layer);

  // Matches the dependency name FieldVerticalInterp declares, and the expression
  // that builds this diag.
  m_diag_name = target.index_name(layer);
  m_mask_name = target.mask_name(layer);

  m_field_in_names.push_back(m_coord_name);
}

void VerticalLevelIndex::
initialize_impl ()
{
  const auto& x = m_fields_in.at(m_coord_name);
  const auto& fid = x.get_header().get_identifier();

  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()==2,
      "Error! Field rank not supported by VerticalLevelIndex.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! VerticalLevelIndex diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // The output is (ncol): the index of the first level past the target in each column
  const auto col_layout = layout.clone().strip_dim(tag);
  FieldIdentifier d_fid (m_diag_name, col_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  m_diagnostic_output = Field(d_fid,true);

  // The mask of columns where the target is valid
  m_diagnostic_output.create_valid_mask(m_mask_name);
}

void VerticalLevelIndex::compute_impl()
{
  using KT = KokkosTypes<DefaultDevice>;

  const Field& x = m_fields_in.at(m_coord_name);
  const auto x_v = x.get_view<const Real**>();
  const auto& xl = x.get_header().get_identifier().get_layout();
  const int ncols = xl.dim(0);
  const int nlevs = xl.dim(1);

  const auto tgt = m_target;
  const auto descending = m_descending;
  const auto extrapolate_bottom = m_extrapolate_bottom;

  auto idx  = m_diagnostic_output.get_view<int*>();
  auto mask = m_diagnostic_output.get_valid_mask().get_view<int*>();

  auto policy = KT::RangePolicy(0,ncols);
  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
    auto beg = ekat::subview(x_v,icol).data();
    auto last = beg + (nlevs-1);

    // beg and last point to the top and bottom entries. The target is before the top entry if
    // it is smaller than it (pressure), or larger than it (heights). It is past the bottom entry
    // if it is larger than it (pressure), or smaller than or equal to it (heights, where the
    // "first level below the target" is nlevs for a target equal to the bottom entry).
    // We don't need to search if the target is out of range: we know the index.
    const bool before_top  = descending ? (tgt>*beg)   : (tgt<*beg);
    const bool past_bottom = descending ? (tgt<=*last) : (tgt>*last);
    if (before_top) {
      idx(icol) = 0;
      mask(icol) = 0;
    } else if (past_bottom) {
      idx(icol) = nlevs;
      mask(icol) = extrapolate_bottom ? 1 : 0;
    } else {
      idx(icol) = first_past(beg,nlevs,tgt,descending);
      mask(icol) = 1;
    }
  });
}

} //namespace scream
