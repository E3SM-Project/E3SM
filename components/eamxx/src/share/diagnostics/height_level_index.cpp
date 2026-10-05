#include "height_level_index.hpp"

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

  m_diag_name = "lev_at_" + z_val + units + "_above_" + surf_ref + "_" + layer;

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

  // The output is (ncol,2): the pair of level indices bracketing the
  // target height in each column (or, if the target is out of range,
  // both indices are set to the bottom level for extrapolation, or to
  // -1 to flag the target as invalid; see class doc above).
  auto idx_layout = layout.clone().strip_dim(tag).append_dim(CMP,2,"bracket");
  FieldIdentifier d_fid (m_diag_name, idx_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  m_diagnostic_output = Field(d_fid,true);
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
  auto idx = m_diagnostic_output.get_view<int**>();
  auto extrapolate_bottom = m_extrapolate_bottom;

  auto policy = KT::RangePolicy(0,ncols);
  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
    auto z_i = ekat::subview(z_v,icol);
    auto beg = z_i.data();
    auto end = beg + nlevs;
    auto it = find_first_smaller_z(beg,end,z_tgt);
    if (it==beg) {
      // Target is above the top entry: always flagged as invalid,
      // regardless of the reference (see class doc).
      idx(icol,0) = -1;
      idx(icol,1) = -1;
    } else if (it==end) {
      // Target is below the bottom entry: extrapolate using it, or
      // flag as invalid, depending on the reference (see class doc).
      idx(icol,0) = extrapolate_bottom ? nlevs-1 : -1;
      idx(icol,1) = extrapolate_bottom ? nlevs-1 : -1;
    } else {
      auto pos = it - beg;
      idx(icol,0) = pos-1;
      idx(icol,1) = pos;
    }
  });
}

} //namespace scream
