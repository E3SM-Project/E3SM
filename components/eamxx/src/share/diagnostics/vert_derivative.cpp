#include "vert_derivative.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {

VertDerivative::
VertDerivative(const ekat::Comm &comm, const ekat::ParameterList &params,
               const std::shared_ptr<const AbstractGrid> &grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  m_coord_name = m_params.get<std::string>("coordinate");

  EKAT_REQUIRE_MSG (m_coord_name!=m_field_name,
      "Error! VertDerivative of a field with respect to itself.\n"
      " - field name: " + m_field_name + "\n");

  m_field_in_names.push_back(m_field_name);
  m_field_in_names.push_back(m_coord_name);
}

void VertDerivative::initialize_impl()
{
  using namespace ShortFieldTagsNames;

  const auto& f      = m_fields_in.at(m_field_name);
  const auto& c      = m_fields_in.at(m_coord_name);
  const auto& fid    = f.get_header().get_identifier();
  const auto& cid    = c.get_header().get_identifier();
  const auto& layout = fid.get_layout();

  // TODO: support higher-dimensioned input fields
  EKAT_REQUIRE_MSG(layout.rank()==2 and (layout.tags().back()==LEV or layout.tags().back()==ILEV),
      "Error! VertDerivative expects a (column,level) field.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(cid.get_layout()==layout,
      "Error! VertDerivative expects the coordinate on the same layout as the field.\n"
      " - field name       : " + fid.name() + "\n"
      " - field layout     : " + layout.to_string() + "\n"
      " - coordinate name  : " + cid.name() + "\n"
      " - coordinate layout: " + cid.get_layout().to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.dims().back()>=2,
      "Error! VertDerivative needs at least two levels.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  auto diag_name = m_field_name + "_differentiate_" + m_coord_name;
  auto d_fid = fid.clone(diag_name).reset_units(fid.get_units()/cid.get_units());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  if (f.has_valid_mask() or c.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }
}

void VertDerivative::compute_impl()
{
  const auto& f = m_fields_in.at(m_field_name);
  const auto& c = m_fields_in.at(m_coord_name);
  const auto& o = m_diagnostic_output;

  // TODO: support higher-dimensioned input fields
  auto f2d = f.get_view<const Real **>();
  auto x2d = c.get_view<const Real **>();
  auto o2d = o.get_view<Real **>();

  using KT          = KokkosTypes<DefaultDevice>;
  using MT          = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;
  using mview_t     = Field::view_dev_t<int**>;
  using cmview_t    = Field::view_dev_t<const int**>;
  const int ncols   = f.get_header().get_identifier().get_layout().dim(0);
  const int nlevs   = f.get_header().get_identifier().get_layout().dim(1);
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  const bool masked   = o.has_valid_mask();
  const bool f_masked = f.has_valid_mask();
  const bool c_masked = c.has_valid_mask();
  auto d_mask = masked   ? o.get_valid_mask().get_view<int**>() : mview_t{};
  auto f_mask = f_masked ? f.get_valid_mask().get_view<const int**>() : cmview_t{};
  auto c_mask = c_masked ? c.get_valid_mask().get_view<const int**>() : cmview_t{};
  const int last = nlevs-1;
  constexpr auto fv = constants::fill_value<Real>;
  auto lambda = KOKKOS_LAMBDA(const MT &team) {
    const int icol = team.league_rank();
    auto fc = ekat::subview(f2d, icol);
    auto xc = ekat::subview(x2d, icol);
    auto oc = ekat::subview(o2d, icol);

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&](const int k) {
      // The levels the stencil reads: k-1,k,k+1 inside, k,k+1 at the top, k-1,k at the bottom
      const int km = k==0    ? 0    : k-1;
      const int kp = k==last ? last : k+1;

      if (masked) {
        int valid = 1;
        for (int j=km; j<=kp; ++j) {
          if (f_masked) valid *= f_mask(icol,j);
          if (c_masked) valid *= c_mask(icol,j);
        }
        d_mask(icol,k) = valid;
        if (valid==0) {
          // Remove this when IO solely relies on mask fields
          oc(k) = fv;
          return;
        }
      }

      if (k==0 or k==last) {
        // One-sided, first order (numpy's gradient with edge_order=1)
        oc(k) = (fc(kp) - fc(km)) / (xc(kp) - xc(km));
      } else {
        // Second order on uneven spacing (numpy's gradient)
        const Real h1 = xc(k) - xc(k-1);
        const Real h2 = xc(k+1) - xc(k);
        oc(k) = ( -h2/(h1*(h1+h2))   * fc(k-1)
                 + (h2-h1)/(h1*h2)   * fc(k)
                 +  h1/(h2*(h1+h2))  * fc(k+1) );
      }
    });
  };
  Kokkos::parallel_for("Compute vertical derivative for " + o.name(), policy, lambda);
}

} // namespace scream
