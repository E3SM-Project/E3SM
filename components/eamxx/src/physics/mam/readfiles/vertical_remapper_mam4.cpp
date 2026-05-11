#include "vertical_remapper_mam4.hpp"
#include <mam4xx/mam4.hpp>

namespace scream
{
VerticalRemapperMAM4::
VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid,
                  const VertRemapType& vremp_type)
:VerticalRemapper(src_grid,
                  tgt_grid) {
  m_vremap_type=vremp_type;
}

void VerticalRemapperMAM4::remap_fwd_impl ()
{
  using namespace ShortFieldTagsNames;

  auto src_vtag = m_src_grid->get_vkind()==AbstractGrid::VKind::Pressure ? LEVP : LEV;

  const auto& src_p = m_src_pressure.at(src_vtag);
  const auto& tgt_p = m_tgt_pressure.at(LEV);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];

    apply_vertical_interpolation(f_src, f_tgt, src_p, tgt_p);
  }
}
/* Invokes MAM4XX routines for vertical interpolation.*/
void VerticalRemapperMAM4::
apply_vertical_interpolation(const Field& f_src, const Field& f_tgt,
                             const Field& p_src, const Field& p_tgt) const
{
  const auto p_tgt_c = p_tgt.get_view<const Real **>();
  const auto datain = f_src.get_view<Real **>();
  const auto dataout =  f_tgt.get_view<Real **>();

  const auto& f_tgt_l = f_tgt.get_header().get_identifier().get_layout();
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = f_tgt_l.dims().back();

  const auto& f_src_l = f_src.get_header().get_identifier().get_layout();
  const int levsiz = f_src_l.dims().back();

  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  const auto policy =
      TPF::get_default_team_policy(ncols, nlevs_tgt);
  using Team = Kokkos::TeamPolicy<KT::ExeSpace>::member_type;
  if (m_vremap_type== MAM4_PSRef) {

    constexpr int unit_factor_pin=1;
    const auto p_src_c = p_src.get_view<const Real **>();

    Kokkos::parallel_for(
      "vert_interp_psref", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        const auto pin_at_icol     = ekat::subview(p_src_c, icol);
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);
        mam4::vertical_interpolation::vert_interp(
            team, levsiz, nlevs_tgt, pin_at_icol, pmid_at_icol, datain_at_icol,
            dataout_at_icol, unit_factor_pin);
          });
  } else if (m_vremap_type== MAM4_ZONAL)
  {
    // unit conversion from mbar->pascals
    const int unit_factor_pin=100;
    const auto p_src_c = p_src.get_view<const Real *>();

    Kokkos::parallel_for(
      "vert_interp_mam4_zonal", policy,
      KOKKOS_LAMBDA(const Team &team) {
        const int icol = team.league_rank();
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);
        mam4::vertical_interpolation::vert_interp(
            team, levsiz, nlevs_tgt, p_src_c, pmid_at_icol, datain_at_icol,
            dataout_at_icol, unit_factor_pin);
          });
  }

}

} // namespace scream
