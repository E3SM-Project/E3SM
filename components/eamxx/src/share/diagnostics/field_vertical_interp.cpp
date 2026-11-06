#include "field_vertical_interp.hpp"
#include "level_index_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>
#include <ekat_units.hpp>

namespace
{
// Linear interpolation (or extrapolation) at c=tgt, of the values f0 and f1, defined at c0 and c1.
// The two forms are algebraically equivalent, but the pressure and the heights have
// historically used a different one, and we keep it to not change answers at roundoff level.
KOKKOS_INLINE_FUNCTION
scream::Real lin_interp (const bool pressure_form,
                         const scream::Real f0, const scream::Real f1,
                         const scream::Real c0, const scream::Real c1,
                         const scream::Real tgt)
{
  return pressure_form
       ? f0 + (f1-f0)/(c1-c0) * (tgt-c0)
       : ( (tgt-c0)*f1 + (c1-tgt)*f0 ) / (c1-c0);
}
} // anonymous namespace

namespace scream
{

FieldVerticalInterp::
FieldVerticalInterp (const ekat::Comm& comm, const ekat::ParameterList& params,
                     const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  m_field_name = m_params.get<std::string>("field_name");

  const level_index_utils::VerticalTarget target (m_params.get<std::string>("vertical_coordinate"),
                                                  m_params.get<std::string>("value"),
                                                  m_params.get<std::string>("units"));
  m_target = target.target;
  m_extrapolate_bottom = target.extrapolate_bottom;
  m_pressure_form = target.coord=="p";

  m_diag_name = m_field_name + target.field_suffix();

  m_field_in_names.push_back(m_field_name);

  // We don't know yet which one (mid/int) we need
  m_coord_name_mid = target.coord_field("mid");
  m_coord_name_int = target.coord_field("int");
  m_field_in_names.push_back(m_coord_name_mid);
  m_field_in_names.push_back(m_coord_name_int);

  // Depend on the (shared) level-index diagnostics, so the search for the level index of
  // this target is only done once, even if multiple fields are requested at this same target.
  // We don't know yet which one (mid/int) we need.
  m_index_name_mid = target.index_name("mid");
  m_index_name_int = target.index_name("int");
  m_field_in_names.push_back(m_index_name_mid);
  m_field_in_names.push_back(m_index_name_int);
}

void FieldVerticalInterp::
initialize_impl ()
{
  const auto& f = m_fields_in.at(m_field_name);
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (f.data_type()==DataType::RealType,
      "Error! FieldVerticalInterp only supports Real data type field.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f.data_type()) + "\n");
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=3,
      "Error! Field rank not supported by FieldVerticalInterp.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"
      "NOTE: if you requested something like 'X.mean('col').interp(..)',\n"
      "      you can avoid this error by requesting 'X.interp(..).mean('col')' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldVerticalInterp diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG (layout.dim(layout.rank()-1)>=2,
      "Error! FieldVerticalInterp needs at least 2 levels to interpolate.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // Figure out which coordinate and level-index fields we need
  m_coord_name = tag==LEV ? m_coord_name_mid : m_coord_name_int;
  m_index_name = tag==LEV ? m_index_name_mid : m_index_name_int;

  // All good, create the diag output
  auto d_fid = fid.clone(m_diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);

  // Add a field representing the mask as extra data to the diagnostic field.
  // If the input is masked, the mask depends on the input's own mask, so we need a mask
  // of our own. Otherwise, the mask only says where the target is valid, so it is the same
  // for any field (and any component of a vector field) at this target and layer:
  // we use the one computed by VerticalLevelIndex, broadcast to the layout of the diagnostic.
  // NOTE: masks with the same name must hold the same values (IO and the horizontal
  //       remapper rely on it); this is true by construction, as the data is shared.
  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
  } else {
    const auto& index_mask = m_fields_in.at(m_index_name).get_valid_mask();
    m_diagnostic_output.set_valid_mask(index_mask.broadcast_to(d_fid.get_layout()));
  }
  m_diagnostic_output.get_header().set_may_be_filled(true);

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

void FieldVerticalInterp::compute_impl()
{
  using RangePolicy = typename KokkosTypes<DefaultDevice>::RangePolicy;
  using cmask2d_t = Field::view_dev_t<const int**>;
  using cmask3d_t = Field::view_dev_t<const int***>;
  using dmask2d_t = Field::view_dev_t<int*>;
  using dmask3d_t = Field::view_dev_t<int**>;

  // The coordinate of the field we interpolate, in 2d (ncol,nlevs)
  const Field& c = m_fields_in.at(m_coord_name);
  const auto c_view = c.get_view<const Real**>();
  const int nlevs = c.get_header().get_identifier().get_layout().dim(1);

  const Field& f = m_fields_in.at(m_field_name);
  const auto& fl = f.get_header().get_identifier().get_layout();

  // The index of the first level past the target in each column, and the mask of columns
  // where the target is valid (both computed once, and shared across all diags targeting
  // this same target, by the VerticalLevelIndex diagnostic).
  // The interpolation uses the levels (k-1,k), with k clamped to [1,nlevs-1]. If k==nlevs
  // in a valid column, we are asked to extrapolate (if the coordinate allows it), and we
  // use the bottom level. Where the mask is 0, the levels are still valid, so we can access
  // data with them, even if we then ignore it.
  const Field& idx_f = m_fields_in.at(m_index_name);
  const auto idx_v = idx_f.get_view<const int*>();
  const auto idx_mask = idx_f.get_valid_mask().get_view<const int*>();

  const auto tgt = m_target;
  const bool extrapolate_bottom = m_extrapolate_bottom;
  const bool pressure_form = m_pressure_form;

  // Whether the INPUT field has its own mask. If not, the mask of our output is the
  // (read-only) shared index mask, which we don't need to compute.
  const bool f_masked = f.has_valid_mask();
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
        const bool extrapolate = extrapolate_bottom and k>=nlevs;
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
            d_view(i) = lin_interp(pressure_form,f_view(i,k0),f_view(i,k1),
                                   c_view(i,k0),c_view(i,k1),tgt);
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
        const bool extrapolate = extrapolate_bottom and k>=nlevs;
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
            d_view(i,j) = lin_interp(pressure_form,f_view(i,j,k0),f_view(i,j,k1),
                                     c_view(i,k0),c_view(i,k1),tgt);
          }
        }
    });
  }

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  auto& mask = m_diagnostic_output.get_valid_mask();
  m_diagnostic_output.deep_copy(constants::fill_value<Real>,mask,true);
}

} //namespace scream
