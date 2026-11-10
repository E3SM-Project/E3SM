#include "pressure_level_index.hpp"
#include "level_index_utils.hpp"

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
  const auto p_value = m_params.get<std::string>("pressure_value");
  const auto units = m_params.get<std::string>("pressure_units");
  m_pressure_level = level_index_utils::pressure_in_Pa(p_value,units);

  const auto& layer = m_params.get<std::string>("vertical_layer");
  EKAT_REQUIRE_MSG (layer=="mid" or layer=="int",
      "Error! Invalid vertical layer for PressureLevelIndex.\n"
      " - input value : " + layer + "\n"
      " - valid values: 'mid', 'int'\n");
  m_pressure_name = "p_" + layer;

  // Matches the dependency name FieldAtPressureLevel declares, and the expression
  // that builds this diag.
  m_diag_name = m_pressure_name + ".searchsorted(" + p_value + ",units='" + units + "')";
  // Named after the request, so that IO names the avg count after the level the user asked for
  m_mask_name = p_value + units + (layer=="int" ? "_int" : "") + "_mask";

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

  // The output is (ncol): the index of the first level with pressure strictly greater
  // than the target, in each column (numpy's searchsorted, with side='right').
  const auto col_layout = layout.clone().strip_dim(tag);
  FieldIdentifier d_fid (m_diag_name, col_layout, ekat::units::none, m_grid->name(), DataType::IntType);
  m_diagnostic_output = Field(d_fid,true);

  // The mask of columns where the target pressure is in range
  m_diagnostic_output.create_valid_mask(m_mask_name);
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

  auto idx  = m_diagnostic_output.get_view<int*>();
  auto mask = m_diagnostic_output.get_valid_mask().get_view<int*>();

  auto policy = KT::RangePolicy(0,ncols);
  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
    auto x1 = ekat::subview(p_v,icol);
    auto beg = x1.data();
    auto end = beg + nlevs;
    auto last = beg + (nlevs-1);
    // The mask says if the target is in range (ends included). The index alone can't
    // tell: nlevs is returned both above the last level and when the target is
    // exactly the last entry.
    // We don't need to search if the target is out of range: we know the index.
    if (p_tgt<*beg) {
      idx(icol) = 0;
      mask(icol) = 0;
    } else if (p_tgt>*last) {
      idx(icol) = nlevs;
      mask(icol) = 0;
    } else {
      // ekat::upper_bound returns the first entry strictly greater than p_tgt
      // (or end, if there is none, which here means p_tgt==*last)
      idx(icol) = ekat::upper_bound(beg,end,p_tgt) - beg;
      mask(icol) = 1;
    }
  });
}

} //namespace scream
