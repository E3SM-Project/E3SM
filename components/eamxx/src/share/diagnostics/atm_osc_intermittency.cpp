// ============================================================
//  atm_osc_intermittency.cpp
// ============================================================
#include "atm_osc_intermittency.hpp"
#include "ekat/ekat_units.hpp"
#include <ekat/ekat_assert.hpp>

namespace scream {

AtmOscIntermittency::AtmOscIntermittency(const ekat::Comm& comm,
                                          const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm, params)
{
  EKAT_REQUIRE_MSG(params.isParameter("tendency_name"),
    "Error! AtmOscIntermittency requires 'tendency_name' in its input parameters.\n");

  m_name          = params.get<std::string>("tendency_name");
  m_noise_floor   = std::stod(params.get<std::string>("noise_floor"));
  m_min_alt_count = std::stoi(params.get<std::string>("min_alt_count"));
  m_diag_name = m_name + "_nf" + params.get<std::string>("noise_floor")
              + "_mac" + params.get<std::string>("min_alt_count")
              + "_atm_osc_intermittency";

  EKAT_REQUIRE_MSG(m_noise_floor >= 0,
    "AtmOscIntermittency: noise_floor must be >= 0");
  EKAT_REQUIRE_MSG(m_min_alt_count >= 1,
    "AtmOscIntermittency: min_alt_count must be >= 1");
}

void AtmOscIntermittency::create_requests()
{
  const auto &gname = m_params.get<std::string>("grid_name");
  add_field<Required>(m_name, gname);
}

void AtmOscIntermittency::initialize_impl(const RunType /* run_type */)
{
  using namespace ekat::units;

  const auto& f_in   = get_field_in(m_name);
  const auto& fid    = f_in.get_header().get_identifier();
  const auto& layout = fid.get_layout();
  const auto& gn     = fid.get_grid_name();

  EKAT_REQUIRE_MSG(
      f_in.data_type() == DataType::RealType,
      "Error! AtmOscIntermittency only supports Real data type field.\n"
      " - field name: " + fid.name() + "\n"
      " - field data type: " + e2str(f_in.data_type()) + "\n");

  // Output is a dimensionless binary flag (0.0 or 1.0) per step.
  // Time-averaging by the output manager converts this to osc_fraction.
  FieldIdentifier d_fid(m_diag_name, layout.clone(), Units::nondimensional(), gn);
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  m_f_pre  = f_in.clone(m_name + "_oi_f_pre");
  m_bt_pre = f_in.clone(m_name + "_oi_bt_pre");
  m_n_alt  = f_in.clone(m_name + "_oi_n_alt");

  m_f_pre.deep_copy(f_in);
  m_bt_pre.deep_copy(0.0);
  m_n_alt.deep_copy(0.0);
}

void AtmOscIntermittency::compute_diagnostic_impl()
{
  using KT = KokkosTypes<DefaultDevice>;
  using RP = typename KT::RangePolicy;

  const auto& f_in  = get_field_in(m_name);
  const int   nelem = f_in.get_header().get_identifier().get_layout().size();

  const Real* f_cur  = f_in.get_internal_view_data<const Real>();
  Real*       f_pre  = m_f_pre.get_internal_view_data<Real>();
  Real*       bt_pre = m_bt_pre.get_internal_view_data<Real>();
  Real*       n_alt  = m_n_alt.get_internal_view_data<Real>();
  Real*       diag   = m_diagnostic_output.get_internal_view_data<Real>();

  const Real noise = m_noise_floor;
  const Real K     = static_cast<Real>(m_min_alt_count);

  Kokkos::parallel_for(m_diag_name, RP(0, nelem), KOKKOS_LAMBDA(int i) {
    const Real bt_cur     = f_cur[i] - f_pre[i];
    const Real abs_bt_cur = bt_cur    < 0 ? -bt_cur    : bt_cur;
    const Real abs_bt_pre = bt_pre[i] < 0 ? -bt_pre[i] : bt_pre[i];

    // Sign-alternation: product is negative AND both BTs exceed noise.
    const bool alternating =
        (bt_cur * bt_pre[i] < Real(0))
     && (abs_bt_cur > noise)
     && (abs_bt_pre > noise);

    // Update consecutive-alternation counter.
    n_alt[i] = alternating ? n_alt[i] + Real(1) : Real(0);

    // Emit binary flag: 1.0 if sustained oscillation detected.
    diag[i] = (n_alt[i] >= K) ? Real(1) : Real(0);

    // Advance carry-forward state.
    bt_pre[i] = bt_cur;
    f_pre[i]  = f_cur[i];
  });
}

} // namespace scream
