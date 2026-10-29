#include "aerocom_weights.hpp"

namespace scream {

AeroComWeights::
AeroComWeights(const ekat::Comm &comm,
               const ekat::ParameterList &params,
               const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_position = m_params.get<std::string>("position");
  m_phase    = m_params.get<std::string>("phase","none");
  m_diag_name = expression(m_position,m_phase);

  m_field_in_names = {"cldfrac_tot", "qc", "qi"};
}

std::string AeroComWeights::
expression (const std::string& position, const std::string& phase)
{
  EKAT_REQUIRE_MSG (position=="top" or position=="bot",
      "Error! Invalid AeroCom position.\n"
      " - position: " + position + "\n"
      " - valid positions: 'top', 'bot'\n");
  EKAT_REQUIRE_MSG (phase=="none" or phase=="liq" or phase=="ice",
      "Error! Invalid AeroCom phase.\n"
      " - phase: " + phase + "\n"
      " - valid phases: 'liq', 'ice', 'none'\n");

  return "cldfrac_tot.aerocom_weights('" + position + "'" +
         (phase=="none" ? std::string("") : ",phase='" + phase + "'") + ")";
}

void AeroComWeights::initialize_impl()
{
  using namespace ShortFieldTagsNames;

  const auto layout = m_grid->get_3d_scalar_layout(LEV);
  for (const auto& [name,f] : m_fields_in) {
    const auto& fid = f.get_header().get_identifier();
    EKAT_REQUIRE_MSG (fid.get_layout().congruent(layout) and f.data_type()==DataType::RealType,
        "Error! AeroComWeights needs Real (COL,LEV) inputs.\n"
        " - field name: " + name + "\n"
        " - field layout: " + fid.get_layout().to_string() + "\n"
        " - field data type: " + e2str(f.data_type()) + "\n");
  }

  FieldIdentifier fid(m_diag_name, layout, ekat::units::none, m_grid->name());
  m_diagnostic_output = Field(fid,true);
}

void AeroComWeights::compute_impl()
{
  /* The probability that layer k is the first cloudy layer seen from above
   * is clr(k-1) - clr(k), where clr(k) is the probability that layers 1..k are
   * all clear. Under the maximum-random overlap assumption, clr is a running
   * product: see get_subcolumn_mask in RRTMGP, which uses equation 14 for its
   * subcolumns; this is equation 13, the column counterpart. The first layer
   * is assumed clear. Seen from below, the same holds with the order of the
   * layers reversed. A layer that is not cloudy has no weight, and leaves clr
   * unchanged.
   */
  const auto cld = m_fields_in.at("cldfrac_tot").get_view<const Real**>();
  const auto qc  = m_fields_in.at("qc").get_view<const Real**>();
  const auto qi  = m_fields_in.at("qi").get_view<const Real**>();
  const auto w   = m_diagnostic_output.get_view<Real**>();

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();

  const bool top = m_position=="top";
  const bool liq = m_phase=="liq";
  const bool ice = m_phase=="ice";
  const Real q_thresh   = q_threshold;
  const Real cld_thresh = cldfrac_tot_threshold;

  // One thread per column: clr is a running product over the levels
  using RangePolicy = Kokkos::RangePolicy<typename KokkosTypes<DefaultDevice>::ExeSpace>;
  Kokkos::parallel_for("AeroComWeights", RangePolicy(0,ncols),
                       KOKKOS_LAMBDA(const int icol) {
    Real clr = 1;
    for (int i=0; i<nlevs-1; ++i) {
      // The layer, and its neighbor on the side we look from
      const int k  = top ? i+1 : nlevs-2-i;
      const int nb = top ? k-1 : k+1;
      const Real q = qc(icol,k) + qi(icol,k);
      Real wk = 0;
      if (q>q_thresh and cld(icol,k)>cld_thresh) {
        const Real clr_k = clr * (1 - ekat::impl::max(cld(icol,nb),cld(icol,k))) /
                                 (1 - ekat::impl::min(cld(icol,nb),1-cld_thresh));
        wk = clr - clr_k;
        if (liq) {
          wk *= qc(icol,k)/q;
        } else if (ice) {
          wk *= 1 - qc(icol,k)/q;
        }
        clr = clr_k;
      }
      w(icol,k) = wk;
    }
    // The first layer seen is assumed clear
    w(icol,top ? 0 : nlevs-1) = 0;
  });
}

}  // namespace scream
