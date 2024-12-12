#include "diagnostics/field_at_height.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/util/ekat_units.hpp"

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
    // if (z>=*mid) {
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

FieldAtHeight::
FieldAtHeight (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
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
}

void FieldAtHeight::
set_grids (const std::shared_ptr<const GridsManager> grids_manager)
{
  const auto& gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_field_name,gname);

  // We don't know yet which one we need
  add_field<Required>(m_z_name+"_mid",gname);
  add_field<Required>(m_z_name+"_int",gname);
}

void FieldAtHeight::
initialize_impl (const RunType /*run_type*/)
{
  const auto& f = get_field_in(m_field_name);
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

  // All good, create the diag output
  FieldIdentifier d_fid (m_diag_name,layout.clone().strip_dim(tag),fid.get_units(),fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src = get_fields_in().front();
  const auto& src_atts = src.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

// =========================================================================================
void FieldAtHeight::compute_diagnostic_impl()
{
  const auto z_view = get_field_in(m_z_name + m_z_suffix).get_view<const Real**>();
  const Field& f = get_field_in(m_field_name);
  const auto& fl = f.get_header().get_identifier().get_layout();

  using RangePolicy = typename KokkosTypes<DefaultDevice>::RangePolicy;

  auto z_tgt = m_z;
  auto nlevs = fl.dims().back();
  if (fl.rank()==2) {
    const auto f_view = f.get_view<const Real**>();
    const auto d_view = m_diagnostic_output.get_view<Real*>();

    RangePolicy policy (0,fl.dims()[0]);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int i) {
        auto f_i = ekat::subview(f_view,i);
        auto z_i = ekat::subview(z_view,i);

        auto beg = z_i.data();
        auto end = beg+nlevs;
        auto it = find_first_smaller_z(beg,end,z_tgt);
        if (it==beg) {
          // We just extapolate with first entry
          d_view(i) = f_i(0);
        } else if (it==end) {
          // We just extapolate with last entry
          d_view(i) = f_i(nlevs-1);
        } else {
          auto pos = it-beg;
          auto z0 = z_i(pos-1);
          auto z1 = z_i(pos);
          auto f0 = f_i(pos-1);
          auto f1 = f_i(pos);

          d_view(i) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
        }
    });
  } else {
    const auto f_view = f.get_view<const Real***>();
    const auto d_view = m_diagnostic_output.get_view<Real**>();

    const auto dim0 = fl.dims()[0];
    const auto dim1 = fl.dims()[1];
    RangePolicy policy (0,dim0*dim1);
    Kokkos::parallel_for(policy,
        KOKKOS_LAMBDA(const int idx) {
        const int i = idx / dim1;
        const int j = idx % dim1;
        auto f_ij = ekat::subview(f_view,i,j);
        auto z_i  = ekat::subview(z_view,i);

        auto beg = z_i.data();
        auto end = beg+nlevs;
        auto it = find_first_smaller_z(beg,end,z_tgt);
        if (it==beg) {
          // We just extapolate with first entry
          d_view(i,j) = f_ij(0);
        } else if (it==end) {
          // We just extapolate with last entry
          d_view(i,j) = f_ij(nlevs-1);
        } else {
          auto pos = it-beg;
          auto z0 = z_i(pos-1);
          auto z1 = z_i(pos);
          auto f0 = f_ij(pos-1);
          auto f1 = f_ij(pos);

          d_view(i,j) = ( (z_tgt-z0)*f1 + (z1-z_tgt)*f0 ) / (z1-z0);
        }
    });
  }
}

} //namespace scream
