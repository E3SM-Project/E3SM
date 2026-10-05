#include "field_at_height.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_units.hpp>

namespace scream
{

FieldAtHeight::
FieldAtHeight (const ekat::Comm& comm, const ekat::ParameterList& params,
               const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  auto surf_ref = m_params.get<std::string>("surface_reference");
  EKAT_REQUIRE_MSG(surf_ref == "sealevel" or surf_ref == "surface",
      "Error! Invalid surface reference for FieldAtHeight.\n"
      " -        field name: " + m_field_name + "\n"
      " - surface reference: " + surf_ref + "\n"
      " -     valid options: sealevel, surface\n");
  m_z_name = (surf_ref == "sealevel") ? "z" : "height";

  const auto units = m_params.get<std::string>("height_units");
  EKAT_REQUIRE_MSG (units=="m",
      "Error! Invalid units for FieldAtHeight.\n"
      " - input units: " + units + "\n"
      " - valid units: m\n");

  auto z_val = m_params.get<std::string>("height_value");
  m_z = std::stod(z_val);
  m_diag_name = m_field_name + "_at_" + z_val + units + "_above_" + surf_ref;

  m_field_in_names.push_back(m_field_name);
  m_field_in_names.push_back(m_z_name+"_mid");
  m_field_in_names.push_back(m_z_name+"_int");

  // Depend on the (shared) bracket-index diagnostics, so the search for
  // the level indices bracketing this height is only done once, even if
  // multiple fields are requested at this same height.
  // We don't know yet which one (mid/int) we need.
  m_index_name_mid = "lev_at_" + z_val + units + "_above_" + surf_ref + "_mid";
  m_index_name_int = "lev_at_" + z_val + units + "_above_" + surf_ref + "_int";
  m_field_in_names.push_back(m_index_name_mid);
  m_field_in_names.push_back(m_index_name_int);
}

void FieldAtHeight::
initialize_impl ()
{
  const auto& f = m_fields_in.at(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (f.data_type()==DataType::RealType,
      "Error! FieldAtHeight only supports Real data type field.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=3,
      "Error! Field rank not supported by FieldAtHeight.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"
      "NOTE: if you requested something like 'field_horiz_avg_at_Y',\n"
      "      you can avoid this error by requesting 'fieldX_at_Y_horiz_avg' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtHeight diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // Figure out the z value
  m_z_suffix = tag==LEV ? "_mid" : "_int";
  m_index_name = tag==LEV ? m_index_name_mid : m_index_name_int;

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);

  // Always create our own valid mask: a target above the model top is
  // always flagged as invalid by the height-index diag we depend on
  // (see HeightLevelIndex), regardless of the surface reference, and a
  // target below the bottom entry is too for "above sealevel" heights;
  // plus the input field itself may carry its own mask.
  m_diagnostic_output.create_valid_mask();
  m_diagnostic_output.get_header().set_may_be_filled(true);

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src = m_fields_in.begin()->second;
  const auto& src_atts = src.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

void FieldAtHeight::compute_impl()
{
  const auto z_view = m_fields_in.at(m_z_name + m_z_suffix).get_view<const Real**>();
  const Field& f = m_fields_in.at(m_field_name);
  const auto& fl = f.get_header().get_identifier().get_layout();

  // The pair of level indices bracketing the target height in each column
  // (computed once, and shared across all diags targeting this same
  // height, by the HeightLevelIndex diagnostic).
  const Field& idx_f = m_fields_in.at(m_index_name);
  const auto idx_v = idx_f.get_view<const int**>();

  using RangePolicy = typename KokkosTypes<DefaultDevice>::RangePolicy;
  using cmask2d_t = Field::view_dev_t<const int**>;
  using cmask3d_t = Field::view_dev_t<const int***>;

  auto z_tgt = m_z;
  constexpr auto fval = constants::fill_value<Real>;
  // Whether the INPUT field has its own mask (our own output mask
  // always exists, see initialize_impl, so no need to track that here).
  bool f_masked = f.has_valid_mask();
  if (fl.rank()==2) {
    const auto f_view = f.get_view<const Real**>();
    const auto d_view = m_diagnostic_output.get_view<Real*>();
    const auto d_mask = m_diagnostic_output.get_valid_mask().get_view<int*>();

    const auto f_mask = f_masked ? f.get_valid_mask().get_view<const int**>() : cmask2d_t{};

    RangePolicy policy (0,fl.dims()[0]);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int i) {
        auto k0 = idx_v(i,0);
        if (k0<0) {
          // Out of range, and not extrapolating (see HeightLevelIndex).
          d_view(i) = fval;
          d_mask(i) = 0;
          return;
        }
        auto k1 = idx_v(i,1);
        if (not f_masked or
            (f_mask(i,k0)!=0 and f_mask(i,k1)!=0)) {
          if (k0==k1) {
            // We just extrapolate with the boundary entry
            d_view(i) = f_view(i,k0);
          } else {
            auto z0 = z_view(i,k0);
            auto z1 = z_view(i,k1);
            auto f0 = f_view(i,k0);
            auto f1 = f_view(i,k1);

            d_view(i) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
          }
          d_mask(i) = 1;
        } else {
          d_mask(i) = 0;
        }
    });
  } else {
    const auto f_view = f.get_view<const Real***>();
    const auto d_view = m_diagnostic_output.get_view<Real**>();
    const auto d_mask = m_diagnostic_output.get_valid_mask().get_view<int**>();

    const auto f_mask = f_masked ? f.get_valid_mask().get_view<const int***>() : cmask3d_t{};

    const auto dim0 = fl.dims()[0];
    const auto dim1 = fl.dims()[1];
    RangePolicy policy (0,dim0*dim1);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int idx2) {
        const int i = idx2 / dim1;
        const int j = idx2 % dim1;
        auto k0 = idx_v(i,0);
        if (k0<0) {
          // Out of range, and not extrapolating (see HeightLevelIndex).
          d_view(i,j) = fval;
          d_mask(i,j) = 0;
          return;
        }
        auto k1 = idx_v(i,1);
        if (not f_masked or
            (f_mask(i,j,k0)!=0 and f_mask(i,j,k1)!=0)) {
          if (k0==k1) {
            // We just extrapolate with the boundary entry
            d_view(i,j) = f_view(i,j,k0);
          } else {
            auto z0 = z_view(i,k0);
            auto z1 = z_view(i,k1);
            auto f0 = f_view(i,j,k0);
            auto f1 = f_view(i,j,k1);

            d_view(i,j) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
          }
          d_mask(i,j) = 1;
        } else {
          d_mask(i,j) = 0;
        }
    });
  }

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  auto& mask = m_diagnostic_output.get_valid_mask();
  m_diagnostic_output.deep_copy(constants::fill_value<Real>,mask,true);
}

} //namespace scream
