#include "field_at_pressure_level.hpp"
#include "level_index_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_team_policy_utils.hpp>
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

  // Figure out the pressure value, and convert to Pa (also validates the units)
  auto p_value = m_params.get<std::string>("pressure_value");
  m_pressure_level = level_index_utils::pressure_in_Pa(p_value,units);

  m_diag_name = m_field_name + "_at_" + p_value + units;

  m_field_in_names.push_back(m_field_name);

  // We don't know yet which one we need
  m_field_in_names.push_back("p_mid");
  m_field_in_names.push_back("p_int");

  // Depend on the (shared) level-index diagnostics, so the binary search
  // for the level index of this pressure level is only done once,
  // even if multiple fields are requested at this same pressure level.
  // We don't know yet which one (mid/int) we need.
  m_index_name_mid = "p_mid.searchsorted(" + p_value + ",units='" + units + "')";
  m_index_name_int = "p_int.searchsorted(" + p_value + ",units='" + units + "')";
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
      "NOTE: if you requested something like 'X.mean('col').interp(..)',\n"
      "      you can avoid this error by requesting 'X.interp(..).mean('col')' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtPressureLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  EKAT_REQUIRE_MSG (layout.dim(layout.rank()-1)>=2,
      "Error! FieldAtPressureLevel needs at least 2 levels to interpolate.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);

  m_pressure_name = tag==LEV ? "p_mid" : "p_int";
  m_index_name = tag==LEV ? m_index_name_mid : m_index_name_int;

  // Add a field representing the mask as extra data to the diagnostic field.
  // If the input is masked, the mask depends on the input's own mask, so we need a mask
  // of our own. Otherwise, the mask only says where the pressure level is in range, so it
  // is the same for any field (and any component of a vector field) at this level and layer:
  // we use the one computed by PressureLevelIndex, broadcast to the layout of the diagnostic.
  // NOTE: masks with the same name must hold the same values (IO and the horizontal
  //       remapper rely on it); this is true by construction, as the data is shared.
  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
  } else {
    // Masks with the same name must hold the same values, hence have the same layout (IO keys the
    // avg count by the mask name, and checks it). So the name of the broadcast mask carries the
    // non-column dims, as for any field that is not (ncol): the (ncol) mask is aliased under that name
    // first, and then broadcast.
    auto index_mask = m_fields_in.at(m_index_name).get_valid_mask();
    const auto& d_layout = d_fid.get_layout();
    if (d_layout.rank()>1) {
      auto mask_name = index_mask.name();
      const std::string end = "_mask";
      if (mask_name.size()>end.size() and mask_name.compare(mask_name.size()-end.size(),end.size(),end)==0) {
        mask_name.erase(mask_name.size()-end.size());
      }
      for (int i=1; i<d_layout.rank(); ++i) {
        mask_name += "_" + d_layout.names()[i] + std::to_string(d_layout.dim(i));
      }
      index_mask = index_mask.alias(mask_name + end);
    }
    m_diagnostic_output.set_valid_mask(index_mask.broadcast_to(d_layout));
  }
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
  using dmask2d_t = Field::view_dev_t<int*>;
  using dmask3d_t = Field::view_dev_t<int**>;

  //This is 2D source pressure
  const Field& p_src = m_fields_in.at(m_pressure_name);
  const auto p_src_v = p_src.get_view<const Real**>();
  const Field& f = m_fields_in.at(m_field_name);

  // The index of the first level with pressure greater than the target in each column,
  // and the mask of columns where the target pressure is in range (both computed
  // once, and shared across all diags targeting this same pressure level, by the
  // PressureLevelIndex diagnostic). The interpolation uses the levels (k-1,k), with
  // k clamped to [1,nlevs-1]: this also covers a target equal to the last entry, and
  // where the mask is 0 the levels are still valid, so we can access data with them,
  // even if we then ignore it.
  const Field& idx_f = m_fields_in.at(m_index_name);
  const auto idx_v = idx_f.get_view<const int*>();
  const auto idx_mask = idx_f.get_valid_mask().get_view<const int*>();

  // The setup for interpolation varies depending on the rank of the input field:
  const int rank = f.rank();

  const auto& pl = p_src.get_header().get_identifier().get_layout();
  const int ncols = pl.dim(0);
  const int nlevs = pl.dim(1);

  auto p_tgt = m_pressure_level;
  bool masked = f.has_valid_mask();
  // NOTE: if the input is not masked, the mask of the diag is the (read-only) shared range
  //       mask, which we don't need to compute. Otherwise, it's our own, and we do.
  if (rank==2) {
    auto policy = KT::RangePolicy(0,ncols);
    auto diag = m_diagnostic_output.get_view<Real*>();
    auto dmask = masked ? m_diagnostic_output.get_valid_mask().get_view<int*>() : dmask2d_t{};
    auto fmask = masked ? f.get_valid_mask().get_view<const int**>() : cmask2d_t{};
    auto f_v  = f.get_view<const Real**>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int icol) {
      const int k  = idx_v(icol);
      const int k1 = k<1 ? 1 : (k>nlevs-1 ? nlevs-1 : k);
      const int k0 = k1-1;
      const bool valid = idx_mask(icol) and (not masked or (fmask(icol,k0)!=0 and fmask(icol,k1)!=0));
      if (masked) {
        dmask(icol) = valid;
      }
      if (valid) {
        // k0 and k1 are always distinct, so this also correctly returns the
        // exact value if p_tgt happens to coincide with the pressure at k0 or k1.
        diag(icol) = f_v(icol,k0) + (f_v(icol,k1)-f_v(icol,k0))/(p_src_v(icol,k1)-p_src_v(icol,k0)) * (p_tgt-p_src_v(icol,k0));
      }
    });
  } else if (rank==3) {
    const int ndims = f.get_header().get_identifier().get_layout().get_vector_dim();
    // Let ekat pick the team size: ndims threads per team may be more than the
    // backend allows (e.g., OpenMP with one thread)
    auto policy = ekat::TeamPolicyFactory<KT::ExeSpace>::get_default_team_policy(ncols,ndims);
    auto diag = m_diagnostic_output.get_view<Real**>();
    auto dmask = masked ? m_diagnostic_output.get_valid_mask().get_view<int**>() : dmask3d_t{};
    auto fmask = masked ? f.get_valid_mask().get_view<const int***>() : cmask3d_t{};
    auto f_v  = f.get_view<const Real***>();
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
      const int k  = idx_v(icol);
      const int k1 = k<1 ? 1 : (k>nlevs-1 ? nlevs-1 : k);
      const int k0 = k1-1;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,ndims),[&](const int idim) {
        const bool valid = idx_mask(icol) and (not masked or (fmask(icol,idim,k0)!=0 and fmask(icol,idim,k1)!=0));
        if (masked) {
          dmask(icol,idim) = valid;
        }
        if (valid) {
          // k0 and k1 are always distinct, so this also correctly returns the
          // exact value if p_tgt happens to coincide with the pressure at k0 or k1.
          diag(icol,idim) = f_v(icol,idim,k0) + (f_v(icol,idim,k1)-f_v(icol,idim,k0))/(p_src_v(icol,k1)-p_src_v(icol,k0)) * (p_tgt-p_src_v(icol,k0));
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
