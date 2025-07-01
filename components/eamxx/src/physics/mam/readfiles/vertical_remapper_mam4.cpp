#include "vertical_remapper_mam4.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include <ekat/util/ekat_units.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>

#include <numeric>

namespace scream
{
VerticalRemapperMAM4::
VerticalRemapperMAM4 (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid,
                  const VertRemapType& vremp_type)
{
  // We only go in one direction for simplicity, since we need to setup some
  // infrsatructures, and we don't want to setup 2x as many "just in case".
  // If you need to remap bwd, just create another remapper with src/tgt grids swapped.
  m_bwd_allowed = false;

  m_vremap_type=vremp_type;

  EKAT_REQUIRE_MSG (src_grid->get_2d_scalar_layout().congruent(tgt_grid->get_2d_scalar_layout()),
      "Error! Source and target grid can only differ for their number of level.\n");

  this->set_grids (src_grid,tgt_grid);
}


void VerticalRemapperMAM4::remap_fwd_impl ()
{
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    apply_vertical_interpolation(f_src, f_tgt, m_src_pmid, m_tgt_pmid);

  }
}
void VerticalRemapperMAM4::
set_source_pressure (const Field& p)
{
  m_src_pmid=p;
}
void VerticalRemapperMAM4::
set_target_pressure (const Field& p)
{
  m_tgt_pmid=p;
}

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

  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncols, nlevs_tgt);

  if (m_vremap_type== MAM4_PSRef) {

    const int unit_factor_pin=1;
    const auto p_src_c = p_src.get_view<const Real **>();

    Kokkos::parallel_for(
      "vert_interp", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
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
      "vert_interp", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol = team.league_rank();
        const auto pmid_at_icol    = ekat::subview(p_tgt_c, icol);
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);
        mam4::vertical_interpolation::vert_interp(
            team, levsiz, nlevs_tgt, p_src_c, pmid_at_icol, datain_at_icol,
            dataout_at_icol, unit_factor_pin);
          });
  } else if (m_vremap_type == MAM4_ELEVATED_EMISSIONS)
  {
      const auto src_x = p_src.get_view<const Real *>();
      const int pverp=nlevs_tgt+1;
      //FIXME: get this values from grid.
      constexpr int nlev = mam4::nlev;
      constexpr Real m2km    = 1e-3;
      Kokkos::parallel_for(
      "tracer_vert_interp_loop", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol = team.league_rank();
        const auto datain_at_icol  = ekat::subview(datain, icol);
        const auto dataout_at_icol = ekat::subview(dataout, icol);

        // FIXME: Try to avoid copy of trg_x by modifying rebin
        // trg_x
        Real trg_x[nlev + 1];
        // I am trying to do this:
        // model_z(1:pverp) = m2km * state(c)%zi(i,pverp:1:-1)
        for(int i = 0; i < pverp; ++i) {
          trg_x[pverp - i - 1] = m2km * p_tgt_c(icol, i);
        }
        team.team_barrier();
        mam4::vertical_interpolation::rebin(team, levsiz, nlevs_tgt,
                                            src_x, trg_x, datain_at_icol, dataout_at_icol);
      });
  }

}

} // namespace scream
