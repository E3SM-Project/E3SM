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

  // Depend on the (shared) level-index diagnostics, so the search for
  // the level index of this height is only done once, even if
  // multiple fields are requested at this same height.
  // We don't know yet which one (mid/int) we need.
  m_index_name_mid = m_z_name + "_mid.searchsorted(" + z_val + ",units='" + units + "')";
  m_index_name_int = m_z_name + "_int.searchsorted(" + z_val + ",units='" + units + "')";
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
      "NOTE: if you requested something like 'X.mean('col').interp(..)',\n"
      "      you can avoid this error by requesting 'X.interp(..).mean('col')' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtHeight diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // Figure out the z value
  m_z_suffix = tag==LEV ? "_mid" : "_int";
  m_index_name = tag==LEV ? m_index_name_mid : m_index_name_int;

  EKAT_REQUIRE_MSG (layout.dim(layout.rank()-1)>=2,
      "Error! FieldAtHeight needs at least 2 levels to interpolate.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);

  // If the input is masked, the mask depends on the input's own mask, so we need a mask
  // of our own. Otherwise, the mask only says where the height is in range (a target above
  // the model top is always flagged as invalid by the HeightLevelIndex diag we depend on,
  // regardless of the surface reference, and a target below the bottom entry is too for
  // "above sealevel" heights), so it is the same for any field (and any component of a vector
  // field) at this height and layer: we use the one computed by HeightLevelIndex, broadcast
  // to the layout of the diagnostic.
  // NOTE: masks with the same name must hold the same values (IO and the horizontal
  //       remapper rely on it); this is true by construction, as the data is shared.
  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
  } else {
    // Masks with the same name must hold the same values, hence have the same layout (IO keys the
    // valid count by the mask name, and checks it). So the name of the broadcast mask carries the
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

  // The index of the first level below the target height in each column, and the mask
  // of columns where the target is valid (both computed once, and shared across all
  // diags targeting this same height, by the HeightLevelIndex diagnostic).
  // The interpolation uses the levels (k-1,k), with k clamped to [1,nlevs-1]. If k==nlevs
  // (valid only for "above surface" heights), we extrapolate with the bottom level.
  // Where the mask is 0, the levels are still valid, so we can access data with them,
  // even if we then ignore it.
  const Field& idx_f = m_fields_in.at(m_index_name);
  const auto idx_v = idx_f.get_view<const int*>();
  const auto idx_mask = idx_f.get_valid_mask().get_view<const int*>();

  using RangePolicy = typename KokkosTypes<DefaultDevice>::RangePolicy;
  using cmask2d_t = Field::view_dev_t<const int**>;
  using cmask3d_t = Field::view_dev_t<const int***>;
  using dmask2d_t = Field::view_dev_t<int*>;
  using dmask3d_t = Field::view_dev_t<int**>;

  auto z_tgt = m_z;
  const int nlevs = z_view.extent(1);
  // Whether the INPUT field has its own mask. If not, the mask of our output is the
  // (read-only) shared range mask, which we don't need to compute.
  bool f_masked = f.has_valid_mask();
  if (fl.rank()==2) {
    const auto f_view = f.get_view<const Real**>();
    const auto d_view = m_diagnostic_output.get_view<Real*>();
    const auto d_mask = f_masked ? m_diagnostic_output.get_valid_mask().get_view<int*>() : dmask2d_t{};

    const auto f_mask = f_masked ? f.get_valid_mask().get_view<const int**>() : cmask2d_t{};

    RangePolicy policy (0,fl.dims()[0]);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int i) {
        const int k  = idx_v(i);
        const int k1 = k<1 ? 1 : (k>nlevs-1 ? nlevs-1 : k);
        const int k0 = k1-1;
        const bool extrapolate = k>=nlevs;
        const bool valid = idx_mask(i) and (not f_masked or
                           (f_mask(i,k1)!=0 and (extrapolate or f_mask(i,k0)!=0)));
        if (f_masked) {
          d_mask(i) = valid;
        }
        if (valid) {
          if (extrapolate) {
            // We just extrapolate with the boundary entry
            d_view(i) = f_view(i,k1);
          } else {
            auto z0 = z_view(i,k0);
            auto z1 = z_view(i,k1);
            auto f0 = f_view(i,k0);
            auto f1 = f_view(i,k1);

            d_view(i) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
          }
        }
    });
  } else {
    const auto f_view = f.get_view<const Real***>();
    const auto d_view = m_diagnostic_output.get_view<Real**>();
    const auto d_mask = f_masked ? m_diagnostic_output.get_valid_mask().get_view<int**>() : dmask3d_t{};

    const auto f_mask = f_masked ? f.get_valid_mask().get_view<const int***>() : cmask3d_t{};

    const auto dim0 = fl.dims()[0];
    const auto dim1 = fl.dims()[1];
    RangePolicy policy (0,dim0*dim1);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int idx2) {
        const int i = idx2 / dim1;
        const int j = idx2 % dim1;
        const int k  = idx_v(i);
        const int k1 = k<1 ? 1 : (k>nlevs-1 ? nlevs-1 : k);
        const int k0 = k1-1;
        const bool extrapolate = k>=nlevs;
        const bool valid = idx_mask(i) and (not f_masked or
                           (f_mask(i,j,k1)!=0 and (extrapolate or f_mask(i,j,k0)!=0)));
        if (f_masked) {
          d_mask(i,j) = valid;
        }
        if (valid) {
          if (extrapolate) {
            // We just extrapolate with the boundary entry
            d_view(i,j) = f_view(i,j,k1);
          } else {
            auto z0 = z_view(i,k0);
            auto z1 = z_view(i,k1);
            auto f0 = f_view(i,j,k0);
            auto f1 = f_view(i,j,k1);

            d_view(i,j) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
          }
        }
    });
  }

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  auto& mask = m_diagnostic_output.get_valid_mask();
  m_diagnostic_output.deep_copy(constants::fill_value<Real>,mask,true);
}

} //namespace scream
