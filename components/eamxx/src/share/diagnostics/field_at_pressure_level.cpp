#include "field_at_pressure_level.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_units.hpp>

namespace scream
{

// =========================================================================================
FieldAtPressureLevel::
FieldAtPressureLevel (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  m_field_name = m_params.get<std::string>("field_name");

  const auto units = m_params.get<std::string>("pressure_units");
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid units for FieldAtPressureLevel.\n"
      " - input units: " + units + "\n"
      " - valid units: 'mb', 'hPa', 'Pa'\n");

  // Figure out the pressure value, and convert to Pa if needed
  auto p_value = m_params.get<std::string>("pressure_value");

  if (units=="mb" || units=="hPa") {
    m_pressure_level = std::stod(p_value)*100;
  } else {
    m_pressure_level = std::stod(p_value);
  }

  m_diag_name = m_field_name + "_at_" + p_value + units;

  m_field_in_names.push_back(m_field_name);

  // We don't know yet which one we need
  m_field_in_names.push_back("p_mid");
  m_field_in_names.push_back("p_int");

  // Depend on the (shared) bracket-index diagnostics, so the binary search
  // for the level indices bracketing this pressure level is only done once,
  // even if multiple fields are requested at this same pressure level.
  // We don't know yet which one (mid/int) we need.
  m_index_name_mid = "lev_at_" + p_value + units + "_mid";
  m_index_name_int = "lev_at_" + p_value + units + "_int";
  m_field_in_names.push_back(m_index_name_mid);
  m_field_in_names.push_back(m_index_name_int);
}

void FieldAtPressureLevel::
initialize_impl ()
{
  const auto& f = m_fields_in.at(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=3,
      "Error! Field rank not supported by FieldAtPressureLevel.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"
      "NOTE: if you requested something like 'field_horiz_avg_at_Y',\n"
      "      you can avoid this error by requesting 'fieldX_at_Y_horiz_avg' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);

  m_pressure_name = tag==LEV ? "p_mid" : "p_int";
  m_index_name = tag==LEV ? m_index_name_mid : m_index_name_int;

  // Add a field representing the mask as extra data to the diagnostic field.
  m_diagnostic_output.create_valid_mask();
  m_diagnostic_output.get_header().set_may_be_filled(true);

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src = m_fields_in.at(m_field_name);
  const auto& src_atts = src.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

void FieldAtPressureLevel::compute_impl()
{
  using KT = KokkosTypes<DefaultDevice>;
  using MemberType = typename KT::MemberType;
  using cmask2d_t = Field::view_dev_t<const int**>;
  using cmask3d_t = Field::view_dev_t<const int***>;

  //This is 2D source pressure
  const Field& p_src = m_fields_in.at(m_pressure_name);
  const auto p_src_v = p_src.get_view<const Real**>();
  const Field& f = m_fields_in.at(m_field_name);

  // The pair of level indices bracketing the target pressure in each
  // column (computed once, and shared across all diags targeting this
  // same pressure level, by the PressureLevelIndex diagnostic). A value
  // of -1 signals that the target pressure is out of range for that column.
  const Field& idx_f = m_fields_in.at(m_index_name);
  const auto idx_v = idx_f.get_view<const int**>();

  // The setup for interpolation varies depending on the rank of the input field:
  const int rank = f.rank();

  const auto& pl = p_src.get_header().get_identifier().get_layout();
  const int ncols = pl.dim(0);

  auto p_tgt = m_pressure_level;
  constexpr auto fval = constants::fill_value<Real>;
  bool masked = f.has_valid_mask();
  if (rank==2) {
    auto policy = KT::RangePolicy(0,ncols);
    auto diag = m_diagnostic_output.get_view<Real*>();
    auto dmask = m_diagnostic_output.get_valid_mask().get_view<int*>();
    auto fmask = masked ? f.get_valid_mask().get_view<const int**>() : cmask2d_t{};
    auto f_v  = f.get_view<const Real**>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
      auto k0 = idx_v(icol,0);
      if (k0<0) {
        diag(icol) = fval;
        dmask(icol) = 0;
        return;
      }
      auto k1 = idx_v(icol,1);
      if (not masked or (fmask(icol,k0)!=0 and fmask(icol,k1)!=0)) {
        // k0 and k1 are always distinct (see PressureLevelIndex), so this
        // also correctly returns the exact value if p_tgt happens to
        // coincide with the pressure at k0 or k1.
        diag(icol) = f_v(icol,k0) + (f_v(icol,k1)-f_v(icol,k0))/(p_src_v(icol,k1)-p_src_v(icol,k0)) * (p_tgt-p_src_v(icol,k0));
        dmask(icol) = 1;
      } else {
        dmask(icol) = 0;
      }
    });
  } else if (rank==3) {
    const int ndims = f.get_header().get_identifier().get_layout().get_vector_dim();
    auto policy = KT::TeamPolicy(ncols,ndims);
    auto diag = m_diagnostic_output.get_view<Real**>();
    auto dmask = m_diagnostic_output.get_valid_mask().get_view<int**>();
    auto fmask = masked ? f.get_valid_mask().get_view<const int***>() : cmask3d_t{};
    auto f_v  = f.get_view<const Real***>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const MemberType& team) {
      int icol = team.league_rank();
      auto k0 = idx_v(icol,0);
      auto k1 = idx_v(icol,1);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims),[&](const int idim) {
        if (k0<0) {
          diag(icol,idim) = fval; // TODO: don't bother setting an arbitrary value
          dmask(icol,idim) = 0;
        } else if (not masked or (fmask(icol,idim,k0)!=0 and fmask(icol,idim,k1)!=0)) {
          // k0 and k1 are always distinct (see PressureLevelIndex), so
          // this also correctly returns the exact value if p_tgt happens
          // to coincide with the pressure at k0 or k1.
          diag(icol,idim) = f_v(icol,idim,k0) + (f_v(icol,idim,k1)-f_v(icol,idim,k0))/(p_src_v(icol,k1)-p_src_v(icol,k0)) * (p_tgt-p_src_v(icol,k0));
          dmask(icol,idim) = 1;
        } else {
          dmask(icol,idim) = 0;
        }
      });
    });
  } else {
    EKAT_ERROR_MSG("Error! field at pressure level only supports fields ranks 2 and 3 \n");
  }

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  auto& mask = m_diagnostic_output.get_valid_mask();
  m_diagnostic_output.deep_copy(constants::fill_value<Real>,mask,true);

}

} //namespace scream
