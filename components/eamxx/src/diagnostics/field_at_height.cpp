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
  int count = end - beg;
  while (count>1) {
    count = end - beg;
    auto mid = beg + count/2;
    if (z>=*mid) {
      beg = mid;
    } else {
      end = mid;
    }
  }
  return z < *beg ? beg : end;
}

} // anonymous namespace

namespace scream
{

// =========================================================================================
FieldAtHeight::
FieldAtHeight (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  const auto& f = params.get<Field>("Field");
  const auto& fid = f.get_header().get_identifier();
  m_field_name = f.name();

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
      " - field layout: " + to_string(layout) + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtHeight diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + to_string(layout) + "\n");

  // Note: you may ask why we can't just store f and be done, rather than go through the
  // add_field infrastructure. Unfortunately, there are some checks in the base classes
  // that require the diagnostic to have 1+ required fields. So we have to do this.
  // TODO: one day we may make atm diags *not* inherit from atm process...
  add_field<Required>(fid);

  // We can also create the diagnostic already!
  const auto& location = params.get<std::string>("Field Level Location");
  const auto diag_field_name = m_field_name + "_at_" + location;

  FieldIdentifier d_fid (diag_field_name,layout.strip_dim(tag),fid.get_units(),fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Figure out the pressure value
  auto chars_start = location.find_first_not_of("0123456789.");
  EKAT_REQUIRE_MSG (chars_start!=0 && chars_start!=std::string::npos,
      "Error! Invalid string for pressure value for FieldAtHeight.\n"
      " - input string   : " + location + "\n"
      " - expected format: Nm, with N integer\n");
  const auto z_str = location.substr(0,chars_start);
  m_z = std::stod(z_str);

  const auto units = location.substr(chars_start);
  EKAT_REQUIRE_MSG (units=="m",
      "Error! Invalid string for height value for FieldAtHeight.\n"
      " - input string   : " + location + "\n"
      " - expected format: Nm, with N integer\n");

  // Create request for z field
  const auto& gname = fid.get_grid_name();
  m_z_name = tag==LEV ? "z_mid" : "z_int";
  add_field<Required>(m_z_name, layout, ekat::units::m, gname);
}

// =========================================================================================
void FieldAtHeight::compute_diagnostic_impl()
{
  const auto z_view = get_field_in(m_z_name).get_view<const Real**>();
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
          // We just extapolate *beg
          d_view(i) = *beg;
        } else if (it==end) {
          // We just extapolate *end
          d_view(i) = *(--end);
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
          d_view(i,j) = *beg;
        } else if (it==end) {
          d_view(i,j) = *(--end);
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
