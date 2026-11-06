#include "height_level_index.hpp"
#include "level_index_utils.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_units.hpp>

namespace
{
// Find first position in array pointed by [beg,end) that is below z
// If all z's in array are >=z, return end
template<typename T>
KOKKOS_INLINE_FUNCTION
const T* find_first_smaller_z (const T* beg, const T* end, const T& z)
{
  // It's easier to find the last entry that is not smaller than z,
  // and then we'll return the ptr after that
  int count = end - beg;
  while (count>1) {
    auto mid = beg + count/2 - 1;
    if (*mid>=z) {
      beg = mid+1;
    } else {
      end = mid+1;
    }
    count = end - beg;
  }

  return *beg < z ? beg : end;
}

} // anonymous namespace

namespace scream
{

HeightLevelIndex::
HeightLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                  const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  auto surf_ref = m_params.get<std::string>("surface_reference");
  EKAT_REQUIRE_MSG(surf_ref == "sealevel" or surf_ref == "surface",
      "Error! Invalid surface reference for HeightLevelIndex.\n"
      " - surface reference: " + surf_ref + "\n"
      " -     valid options: sealevel, surface\n");
  const std::string z_name = (surf_ref == "sealevel") ? "z" : "height";
  // An "above surface" height is always well defined near the surface
  // (every column has one), so extrapolate a target below the bottom
  // entry. An "above sealevel" height is not (e.g. a target elevation
  // below a mountain's surface), so flag it as invalid instead. Either
  // way, a target above the top entry is always flagged as invalid; see
  // class doc.
  m_extrapolate_bottom = (surf_ref == "surface");

  const auto units = m_params.get<std::string>("height_units");
  EKAT_REQUIRE_MSG (units=="m",
      "Error! Invalid units for HeightLevelIndex.\n"
      " - input units: " + units + "\n"
      " - valid units: m\n");

  auto z_val = m_params.get<std::string>("height_value");
  m_z = std::stod(z_val);

  const auto& layer = m_params.get<std::string>("vertical_layer");
  EKAT_REQUIRE_MSG (layer=="mid" or layer=="int",
      "Error! Invalid vertical layer for HeightLevelIndex.\n"
      " - input value : " + layer + "\n"
      " - valid values: 'mid', 'int'\n");
  m_z_field_name = z_name + "_" + layer;

  // Matches the dependency name FieldAtHeight declares, and the expression
  // that builds this diag.
  m_diag_name = m_z_field_name + ".searchsorted(" + z_val + ",units='" + units + "')";
  // Named after the request, so that IO names the avg count after the height the user asked for
  m_mask_name = z_val + units + "_above_" + surf_ref + (layer=="int" ? "_int" : "") + "_mask";

  m_field_in_names.push_back(m_z_field_name);
}

void HeightLevelIndex::
initialize_impl ()
{
  const auto& z = m_fields_in.at(m_z_field_name);
  const auto& fid = z.get_header().get_identifier();

  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()==2,
      "Error! Field rank not supported by HeightLevelIndex.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! HeightLevelIndex diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // The output is (ncol): the index of the first level with height strictly below the target
  // in each column (nlevs if there is none). Since heights decrease with the level index,
  // this is searchsorted's counterpart (with side='right') for a descending array.
  const auto col_layout = layout.clone().strip_dim(tag);
  FieldIdentifier d_fid (m_diag_name, col_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  m_diagnostic_output = Field(d_fid,true);

  // The mask of columns where the target height is valid
  m_diagnostic_output.create_valid_mask(m_mask_name);
}

void HeightLevelIndex::compute_impl()
{
  using KT = KokkosTypes<DefaultDevice>;

  const Field& z = m_fields_in.at(m_z_field_name);
  const auto z_v = z.get_view<const Real**>();
  const auto& zl = z.get_header().get_identifier().get_layout();
  const int ncols = zl.dim(0);
  const int nlevs = zl.dim(1);

  auto z_tgt = m_z;
  auto idx = m_diagnostic_output.get_view<int*>();
  auto mask = m_diagnostic_output.get_valid_mask().get_view<int*>();
  auto extrapolate_bottom = m_extrapolate_bottom;

  auto policy = KT::RangePolicy(0,ncols);
  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
    auto z_i = ekat::subview(z_v,icol);
    auto beg = z_i.data();
    auto end = beg + nlevs;
    auto last = beg + (nlevs-1);
    // Heights decrease with the level index, so *beg is the top, and *last the bottom.
    // We don't need to search if the target is out of range: we know the index.
    if (z_tgt>*beg) {
      // Target is above the top entry: always invalid, regardless of the reference.
      idx(icol) = 0;
      mask(icol) = 0;
    } else if (z_tgt<=*last) {
      // Target is at or below the bottom entry: extrapolate, or flag it as invalid,
      // depending on the reference. See class doc.
      idx(icol) = nlevs;
      mask(icol) = extrapolate_bottom ? 1 : 0;
    } else {
      idx(icol) = find_first_smaller_z(beg,end,z_tgt) - beg;
      mask(icol) = 1;
    }
  });
}

} //namespace scream
