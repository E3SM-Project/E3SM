#include "pressure_level_index.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_upper_bound.hpp>
#include <ekat_units.hpp>

namespace scream
{

// =========================================================================================
PressureLevelIndex::
PressureLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  const auto units = m_params.get<std::string>("pressure_units");
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid units for PressureLevelIndex.\n"
      " - input units: " + units + "\n"
      " - valid units: 'mb', 'hPa', 'Pa'\n");

  // Figure out the pressure value, and convert to Pa if needed
  auto p_value = m_params.get<std::string>("pressure_value");

  if (units=="mb" || units=="hPa") {
    m_pressure_level = std::stod(p_value)*100;
  } else {
    m_pressure_level = std::stod(p_value);
  }

  const auto& layer = m_params.get<std::string>("vertical_layer");
  EKAT_REQUIRE_MSG (layer=="mid" or layer=="int",
      "Error! Invalid vertical layer for PressureLevelIndex.\n"
      " - input value : " + layer + "\n"
      " - valid values: 'mid', 'int'\n");
  m_pressure_name = "p_" + layer;

  m_diag_name = "lev_at_" + p_value + units + "_" + layer;

  m_field_in_names.push_back(m_pressure_name);
}

void PressureLevelIndex::
initialize_impl ()
{
  const auto& p = m_fields_in.at(m_pressure_name);
  const auto& fid = p.get_header().get_identifier();

  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()==2,
      "Error! Field rank not supported by PressureLevelIndex.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! PressureLevelIndex diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // The output is (ncol,2): the pair of level indices bracketing the
  // target pressure in each column.
  auto idx_layout = layout.clone().strip_dim(tag).append_dim(CMP,2,"bracket");
  FieldIdentifier d_fid (m_diag_name, idx_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  m_diagnostic_output = Field(d_fid,true);
}

void PressureLevelIndex::compute_impl()
{
  using KT = KokkosTypes<DefaultDevice>;

  const Field& p = m_fields_in.at(m_pressure_name);
  const auto p_v = p.get_view<const Real**>();
  const auto& pl = p.get_header().get_identifier().get_layout();
  const int ncols = pl.dim(0);
  const int nlevs = pl.dim(1);

  auto p_tgt = m_pressure_level;

  auto idx  = m_diagnostic_output.get_view<int**>();

  auto policy = KT::RangePolicy(0,ncols);
  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
    auto x1 = ekat::subview(p_v,icol);
    auto beg = x1.data();
    auto end = beg + nlevs;
    auto last = beg + (nlevs-1);
    if (p_tgt<*beg or p_tgt>*last) {
      // Target pressure out of range: signal it with a -1 sentinel.
      idx(icol,0) = -1;
      idx(icol,1) = -1;
    } else {
      // ekat::upper_bound returns the index of the first entry strictly
      // greater than p_tgt. Clamp it to [1,nlevs-1], so the bracket
      // [k1-1,k1] is always made of two distinct levels. This also
      // correctly handles p_tgt exactly matching a level's pressure
      // (either endpoint): the general linear-interpolation formula
      // used downstream evaluates exactly to that level's value in
      // that case, so no separate corner case is needed.
      auto ub = ekat::upper_bound(beg,end,p_tgt);
      auto k1 = ub - beg;
      if (k1<1) {
        k1 = 1;
      } else if (k1>nlevs-1) {
        k1 = nlevs-1;
      }
      idx(icol,0) = k1-1;
      idx(icol,1) = k1;
    }
  });
}

} //namespace scream
